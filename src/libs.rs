use std::collections::HashMap;
use std::path::Path;
use rust_htslib::{bam, bam::Read};
use std::str::from_utf8;
use bio::io::fasta::IndexedReader;
use chrono::{DateTime, Local};
use std::error;
use itertools::Itertools;
extern crate bambam;
use std::fs;
use crossbeam::channel::{unbounded};
use log::debug;




// Change the alias to `Box<error::Error>`.
type Result<T> = std::result::Result<T, Box<dyn error::Error>>;






#[derive(Debug,Default,Clone)]
/// structure which holds our
/// counts of nucleotides for 
/// any given genomic position
pub struct Pileup {
    /// chromosome
    chromosome: String,
    /// start is 0 based similar to BED
    start: u64,
	/// end is 0 based
	end: u64,
    /// the reference nuc if 
	/// available through provided 
	/// reference FASTA
    reference: String ,
    /// number of observe adenines
    nuc_a: u32,
    /// number of observe thyosins
    nuc_t: u32,
    /// number of observe cytosins
    nuc_c: u32,
    /// number of observe guanins
    nuc_g: u32,
    /// number of observe ambigious
	/// elements, e.g. N or X
    ambigious: u32,
	/// number of reads with insertion
	ins: u32,
	/// number of reads with deletion
	del: u32,
	/// depth at that given position
	depth: u32,
	/// mutated or not
	mutated: bool,
	/// sum of variant frequency (sum/depth)
	vaf: f32,
	//// reference allele frequency (ref/depth)
	raf: f32,
}

// takes the generates pileups in our specific
// format and prints them according to our wishes
// currently only tsv, might add more possibilities
// in the future
pub fn print_pileup(
	result:&[Pileup],
	out: &str, 
	version: &str, 
	author: &str , 
	command: &str 
)-> Result<i32> {
	let now: DateTime<Local> = Local::now();
	let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(&out)?;
	writer.write_record(&["##","abc:",version,"","","","","","","","","","",""])?;
    writer.write_record(&["##","author:",author,"","","","","","","","","","",""])?;
    writer.write_record(&["##","date:",&now.to_rfc2822(),"","","","","","","","","","",""])?;
    writer.write_record(&["##","command:", command,"","","","","","","","","","",""])?;    
    writer.write_record(&["# chromosome","start","end","reference","A","T","C","G","ambigious","ins","del","depth","VAF","RAF"])?;
    for element in result.iter() {
        writer.write_record(&[
			element.chromosome.clone(),
			element.start.to_string(),
			element.end.to_string(),
			element.reference.clone(),
			element.nuc_a.to_string(),
			element.nuc_t.to_string(),
			element.nuc_c.to_string(),
			element.nuc_g.to_string(),
			element.ambigious.to_string(),
			element.ins.to_string(),
			element.del.to_string(),
			element.depth.to_string(),
			format!("{:.4}",element.vaf),
			format!("{:.4}",element.raf),
			])?;
    };
    writer.flush()?;
    Ok(0)
}


// this simply takes a BED file and organizes entries in a Hashmap
// where they are organized by chromosome which should speedup
// later the indexed access
pub fn parse_bed_file (
	input: &str
) -> HashMap<String,Vec<u64>>{
	assert!(
		Path::new(input).exists(),
		"ERROR: BED file {} does not exist!",
		input,
	);
	let mut bed_result : HashMap<String,Vec<u64>> = HashMap::new();
	let mut bed = bio::io::bed::Reader::from_file(input).unwrap();
	for entry in bed.records(){
		let record = entry.expect("ERROR: wrong BED record");
        // dirty hack, if that happens we create for each position a single entry
        // we should still warn though that this is currently not recommended
		if record.end()-record.start()>1 {
			eprintln!("WARNING: entry {:?} contained not a position but a range! This will cost you.....",&record);
            for pos in record.start()..record.end() {
                bed_result
                .entry(record.chrom().to_string())
                .or_insert_with(Vec::new)
                .push(pos);        
            }
		}
		bed_result
			.entry(record.chrom().to_string())
			.or_insert_with(Vec::new)
			.push(record.start());
	
	};
	for (_, entries) in bed_result.iter_mut(){
		entries.sort_unstable();
		entries.dedup();
	}
	bed_result
}

// This reads now the BAM file and uses the positions
// to check for each position the observed number of nucleotides
// It sorts as well again by chromosome and position to be safe
pub fn analyze_bam_positions (
	input: &str , 
	positions: &HashMap<String,Vec<u64>>, 
	ref_fasta: &str, 
	threads: &usize
) -> Vec<Pileup> {
	assert!(
		Path::new(input).exists(),
		"ERROR: BAM file {} does not exist!",
		input,
	);
	// for the threads we need to have to first clone everything
	// otherwise the "scope" part cant share the common variables
	let input_local     = &input.to_owned();
	let positions_local = positions.clone();
	// less ideal but we have thanks to clap a reference 
	// but want to avoid double referencing, therefore the
	// work-around here
	let ref_local       = <&str>::clone(&ref_fasta);
	//let cpus = 2;
	//let n_chrom =positions_local.len();
	//let chunk_len = (n_chrom as u32/threads) as usize + 1 ;
	//eprintln!("INFO: Foundchromosomes: {} divided by threads {}", n_chrom, chunk_len);
	//let reference_local = ref_fasta.to_owned();
	let thread_results : Vec<Pileup> ;
	let (snd, rxv) = unbounded();
	rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();
	eprintln!("Current thread-pool number: {:?}",rayon::current_num_threads());
	// now we launch for each position a new
	// pileup call and send the result to the receiver
	rayon::scope( |child| {
		for chr in positions_local.keys().sorted() {
			for pos in positions_local.get(chr).unwrap().iter().sorted(){
				// cloning the sender will still send everything into same receiver
				let snd_local = snd.clone();
				// now we define how many threads should be used
				child.spawn( move |_| {
					//eprintln!("Current thread index: {:?}",rayon::current_thread_index());
					let mut bam_file = bam::IndexedReader::from_path(input_local).unwrap();
					let tmp_result = match ref_local { 
						"NONE" => {
									fetch_position(&mut bam_file,chr,pos,None,false)
								}
								,
						_     => {
									let faidx = IndexedReader::from_file(&ref_local).unwrap();
									fetch_position(&mut bam_file,chr,pos,Some(faidx),false)
								}
					};
					snd_local.send(tmp_result).expect("ERROR: thread could not communicate result!");
				});	
			}
		}
	
	});
	// we need to close afterwards 
	// the original sender as it otherwise still waits
	drop(snd);
	thread_results = rxv.iter().collect();
	thread_results
}

// this sub-function gets the iterator result from the positions
// and return the result of pileup analysis at that given position
// allows to potentially already parallelize over each position analyzed
fn fetch_position(
	bam: & mut bam::IndexedReader, 
	chr: &str, 
	pos:&u64, 
	ref_file: Option<IndexedReader<fs::File>>,
	_legacy: bool,
) -> Pileup {
	// NOTE:
	// SAM is 1 based and not 0 based, need to correct for that in
	let start = *pos ;
	let end   = *pos +1 ;
	let mut with_ref = false;
	// this obtains now the pileup at that
	// given position
	bam.fetch((chr,start,end)).expect("ERROR: could not fetch region");
	// currently we ignore completely clipping
	let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
	match ref_file {
		Some(mut x) => {
				with_ref = true;
				x.fetch(chr,*pos,pos+1).expect("ERROR: could not fetch interval on reference ");
				let mut ref_seq = Vec::new();
				x.read(&mut ref_seq).expect("ERROR: could not read sequence from reference");
				collection.reference = String::from(from_utf8(&ref_seq).unwrap())
			}
		None => {collection.reference = String::from("NA"); collection.mutated = false }
	}
	collection.chromosome = chr.to_string();
	collection.start   = start;
	collection.end     = end;
	for pile in bam.pileup().map(|x| x.expect("ERROR: could not parse BAM file")){
		// now we only care about the position we inquire
		//dbg!(&pile);
		if pile.pos() as u64 == start {
			collection.depth = pile.depth();
			for alignment in pile.alignments(){
				// some have none, no idea why.
				// additionally the is Del from indel is not working as anticipated
				// the is_del works instead though...
				if alignment.is_del() {
					collection.del += 1;
				}else{
					match alignment.indel() {
						bam::pileup::Indel::Ins(_len) => {collection.ins += 1; continue } ,
						bam::pileup::Indel::Del(_len) => {collection.del += 1; continue },
						bam::pileup::Indel::None => ()
					}
				}
				
				// sometimes we get reads with 0 length, no idea why
				if alignment.record().seq_len() == 0 { continue }
				let qpos = match  alignment.qpos(){
					Some(q) => q,
					_ =>  continue,
				};
				let nuc  = &alignment.record().seq().as_bytes()[qpos..qpos+1];
				let nuc1 = from_utf8(nuc).unwrap();
				match nuc1.to_uppercase().as_str() {
					"A" => collection.nuc_a += 1  ,
					"T" => collection.nuc_t += 1  ,
					"C" => collection.nuc_c += 1  ,
					"G" => collection.nuc_g += 1  ,
					_ => collection.ambigious += 1  ,
				}
			}
			break;
		}
	};
	if with_ref {
		// if legacy {
		// 	eval_mutation_legacy(collection)
		// }else{
			eval_mutation(collection)
		//}
	}else{
		collection
	}
}


/// this function will take a pileup information in our native
/// format and evaluate the necessary seeked information.
/// It then returns the pileup information with additional
/// filled fields
/// Steps: 
/// - build hash with observations
/// - identify reference nucleotide => RAF
/// - remove element which is ref
/// - identify element with most counts
/// - calculate VAF
fn eval_mutation(
	mut obs:Pileup
) -> Pileup{
	debug!("{:?}",obs);
	let mut dominant_list : HashMap<&str,u32> =
	[
		 ("A", obs.nuc_a), 
		 ("T", obs.nuc_t),
		 ("C", obs.nuc_c),
		 ("G", obs.nuc_g), 
		 ("N", obs.ambigious),
		 ("I", obs.ins),
		 ("D", obs.del) 
		 ]
		.iter()
		.cloned()
		.collect();
	// check that the depth equals sum of observations
	if obs.nuc_a + obs.nuc_t + obs.nuc_c + obs.nuc_g + obs.ambigious + obs.ins + obs.del != obs.depth {
		debug!("A {} T {} C{} G {} N {} INS{} DEL {} Depth {}",obs.nuc_a,obs.nuc_t,obs.nuc_c,obs.nuc_g,obs.ambigious,obs.ins,obs.del,obs.depth );
		panic!("ERROR: the sum of observations does not equal depth !");
	}
	match obs.reference.to_uppercase().as_str(){
		"A" => {
			if obs.nuc_a != obs.depth {
				// get ref nuc
				let original = *dominant_list.get("A").expect("ERROR: A did not exist!!") as f32;
				// remove ref nuc
				dominant_list.remove("A");
				let depth    = (obs.depth) as f32;
				// calc raf
				obs.raf      = original/depth;
				// next we want to get most dominant alternative mutated nucleotid
				let max_nuc  = dominant_list
					.iter()
					.max_by_key(|entry| entry.1).expect("ERROR: no element in nuc hash!");
				// calc vaf of dominant second (if tied, still correct)
				obs.vaf      = *max_nuc.1 as f32/depth;
				// define as being non-reference observation
				obs.mutated  = true;
			}
		},
		"T" => {
			if obs.nuc_t != obs.depth {
				let original = *dominant_list.get("T").expect("ERROR: A did not exist!!") as f32;
				dominant_list.remove("T");
				let depth    = (obs.depth) as f32;
				obs.raf      = original/depth;
				// next we want to get most dominant alternative mutated nucleotide
				let max_nuc  = dominant_list
					.iter()
					.max_by_key(|entry| entry.1).expect("ERROR: no element in nuc hash!");
				obs.vaf      = *max_nuc.1 as f32/depth;
				obs.mutated  = true;
			}
		},
		"C" => {
			if obs.nuc_c != obs.depth {
				let original = *dominant_list.get("C").expect("ERROR: A did not exist!!") as f32;
				dominant_list.remove("C");
				let depth    = (obs.depth) as f32;
				obs.raf      = original/depth;
				// next we want to get most dominant alternative mutated nucleotide
				let max_nuc  = dominant_list
					.iter()
					.max_by_key(|entry| entry.1).expect("ERROR: no element in nuc hash!");
				obs.vaf      = *max_nuc.1 as f32/depth;
				obs.mutated  = true;
			}
		},
		"G" => {
			if obs.nuc_g != obs.depth {
				let original = *dominant_list.get("G").expect("ERROR: A did not exist!!") as f32;
				dominant_list.remove("G");
				let depth    = (obs.depth) as f32;
				obs.raf      = original/depth;
				// next we want to get most dominant alternative mutated nucleotide
				let max_nuc  = dominant_list
					.iter()
					.max_by_key(|entry| entry.1).expect("ERROR: no element in nuc hash!");
				obs.vaf      = *max_nuc.1 as f32/depth;
				obs.mutated  = true;
			}else{
				obs.mutated  = false;
				obs.vaf      = 0.0_f32;
				obs.raf      = 1.0_f32;
			}
		},
		"N" => (),
		_ => {panic!("ERROR: non compatible base found in reference sequence: {}!",obs.reference.as_str())},
	}	 
	
	obs
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
	
    #[test]
	/// test for a 100% reference observations
    fn test_eva_mut_1() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_a = 10;
		collection.depth = 10;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,false);
		assert_eq!(result.raf,1.0_f32);
		assert_eq!(result.vaf,0.0_f32);

    }

	#[test]
	/// test for a 100% mutation observations
    fn test_eva_mut_2() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_t = 10;
		collection.depth = 10;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,true);
		assert_eq!(result.raf,0.0_f32);
		assert_eq!(result.vaf,1.0_f32);

    }


	#[test]
	/// test for a 50% mutation observations
    fn test_eva_mut_3() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_t = 10;
		collection.nuc_a = 10;
		collection.depth = 20;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,true);
		assert_eq!(result.raf,0.5_f32);
		assert_eq!(result.vaf,0.5_f32);

    }


	#[test]
	/// test for a 25% multi allelic mutation observations, 2 with equal proportion
    fn test_eva_mut_4() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_t = 5;
		collection.nuc_c = 5;
		collection.nuc_a = 10;
		collection.depth = 20;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,true);
		assert_eq!(result.raf,0.5_f32);
		assert_eq!(result.vaf,0.25_f32);

    }


	#[test]
	/// test for a 25% multi allelic mutation observations, 3 with unequal proportion
    fn test_eva_mut_5() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_t = 2;
		collection.nuc_g = 3;
		collection.nuc_c = 5;
		collection.nuc_a = 10;
		collection.depth = 20;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,true);
		assert_eq!(result.raf,0.5_f32);
		assert_eq!(result.vaf,0.25_f32);

    }

	#[test]
	#[should_panic(expected = "ERROR: the sum of observations does not equal depth !")]
	/// test that fails if depth is unequal to depth
    fn test_eva_mut_6() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_a = 10;
		collection.depth = 20;
		eval_mutation(collection);
    }

	#[test]
	/// test if indels are present 
    fn test_eva_mut_7() {
        let mut collection : Pileup = Pileup { vaf: 0_f32, raf: 1_f32, ..Default::default() };
		collection.reference = String::from("A");
		collection.nuc_a = 10;
		collection.depth = 20;
		collection.del   = 10;
		let result = eval_mutation(collection);
		assert_eq!(result.mutated,true);
		assert_eq!(result.raf,0.5_f32);
		assert_eq!(result.vaf,0.5_f32);
    }

   
}