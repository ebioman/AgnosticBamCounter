use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::collections::HashMap;
use std::path::Path;
use rust_htslib::{bam, bam::Read};
use std::str::from_utf8;
use bio::io::fasta::IndexedReader;
use chrono::{DateTime, Local};
use std::error;
use std::process;
use std::env;
use itertools::Itertools;
extern crate bambam;
use std::fs;
use crossbeam::channel::{unbounded};
use rayon;





// Change the alias to `Box<error::Error>`.
type Result<T> = std::result::Result<T, Box<dyn error::Error>>;




#[derive(Debug,Default,Clone)]
/// structure which holds our
/// counts of nucleotides for 
/// any given genomic position
struct Pileup {
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
    nuc_a: i32,
    /// number of observe thyosins
    nuc_t: i32,
    /// number of observe cytosins
    nuc_c: i32,
    /// number of observe guanins
    nuc_g: i32,
    /// number of observe ambigious
	/// elements, e.g. N or X
    ambigious: i32,
	/// number of reads with insertion
	ins: i32,
	/// number of reads with deletion
	del: i32,
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
fn print_pileup(result:&[Pileup],out: &str, version: &str, author: &str , command: &str )-> Result<i32> {
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
fn parse_bed_file (input: &str) -> HashMap<String,Vec<u64>>{
	assert!(
		Path::new(input).exists(),
		"ERROR: BED file {} does not exist!",
		input,
	);
	let mut bed_result : HashMap<String,Vec<u64>> = HashMap::new();
	let mut bed = bio::io::bed::Reader::from_file(input).unwrap();
	for entry in bed.records(){
		let record = entry.expect("ERROR: wrong BED record");
		if record.end()-record.start()>1 {
			panic!("ERROR: entry {:?} contained not a position but a range!",&record);
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
fn analyze_bam_positions (input: &str , positions: &HashMap<String,Vec<u64>>, ref_fasta: &str, threads: &usize) -> Vec<Pileup>{
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
									fetch_position(&mut bam_file,chr,pos,None)
								}
								,
						_     => {
									let faidx = IndexedReader::from_file(&ref_local).unwrap();
									fetch_position(&mut bam_file,chr,pos,Some(faidx))
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
fn fetch_position(bam: & mut bam::IndexedReader, chr: &str, pos:&u64, ref_file: Option<IndexedReader<fs::File>>) -> Pileup {
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
				// some have none, no idea why
				match alignment.indel() {
					bam::pileup::Indel::Ins(_len) => {collection.ins += 1; continue } ,
    				bam::pileup::Indel::Del(_len) => {collection.del += 1; continue },
    				bam::pileup::Indel::None => ()
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
		eval_mutation(collection)
	}else{
		collection
	}
}

fn eval_mutation(mut obs:Pileup) -> Pileup{
	match obs.reference.to_uppercase().as_str(){
		"A" => {
			if ( obs.nuc_c !=0 ) || (obs.nuc_g != 0) || (obs.nuc_t != 0) || (obs.del !=0) || (obs.ins !=0) {
				let mutated = (obs.nuc_c + obs.nuc_g + obs.nuc_t  + obs.del + obs.ins) as f32;
				let original = (obs.nuc_a) as f32;
				let depth = (obs.depth) as f32;
				obs.vaf = mutated/depth; 
				obs.raf = original/depth;
				obs.mutated = true;
			}
		},
		"T" => {
			if ( obs.nuc_c !=0 ) || (obs.nuc_g != 0) || (obs.nuc_a != 0) || (obs.del !=0) || (obs.ins !=0){
				let mutated = (obs.nuc_c + obs.nuc_g + obs.nuc_a  + obs.del + obs.ins) as f32;
				let original = (obs.nuc_t) as f32;
				let depth = (obs.depth) as f32;
				obs.vaf = mutated/depth; 
				obs.raf = original/depth;
				obs.mutated = true;
			}
		},
		"C" => {
			if ( obs.nuc_t !=0 ) || (obs.nuc_g != 0) || (obs.nuc_a != 0) || (obs.del !=0) || (obs.ins !=0){
				let mutated = (obs.nuc_a + obs.nuc_g + obs.nuc_t  + obs.del + obs.ins) as f32;
				let original = (obs.nuc_c) as f32;
				let depth = (obs.depth) as f32;
				obs.vaf = mutated/depth; 
				obs.raf = original/depth;
				obs.mutated = true;
			}
		},
		"G" => {
			if ( obs.nuc_t !=0 ) || (obs.nuc_c != 0) || (obs.nuc_a != 0) || (obs.del !=0) || (obs.ins !=0){
				let mutated = (obs.nuc_c + obs.nuc_a + obs.nuc_t  + obs.del + obs.ins) as f32;
				let original = (obs.nuc_g) as f32;
				let depth = (obs.depth) as f32;
				obs.vaf = mutated/depth; 
				obs.raf = original/depth;
				obs.mutated = true;
			}
		},
		"N" => (),
		_ => {panic!("ERROR: non compatible base found in reference sequence: {}!",obs.reference.as_str())},
	}	
	obs
}

fn main() {
	// now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
		.arg(Arg::with_name("BAM")
			.short("b")
			.long("bam")
			.value_name("BAM")
			.help("aligned short or long reads")
			.takes_value(true)
			.required(true))
		.arg(Arg::with_name("POSITIONS")
			.short("p")
			.long("positions")
			.value_name("BED")
			.help("A bed file with positions")
			.takes_value(true)
			.required(true))
        .arg(Arg::with_name("OUT")
			.short("o")
			.long("outfile")
			.value_name("TSV")
			.help("the output file which will be in tsv format")
			.takes_value(true)
			.required(true))
		.arg(Arg::with_name("REF")
			.short("r")
			.long("reference")
			.value_name("FASTA")
			.help("if reference in fasta format provided, reference nucleotide \
                 \nis provided for each positions, otherwise NA")
			.takes_value(true)
            .required(false))
		.arg(Arg::with_name("THREADS")
			.short("t")
			.long("threads")
			.value_name("INT")
			.help("number of threads used in thread-pool for querying positions ")
			.default_value("1")
			.takes_value(true)
			.required(false))
		.arg(Arg::with_name("BAMBAM")
			.short("x")
			.long("bambam")
			.value_name("IDH")
			.takes_value(false)
			.required(false)
			.hidden(true))	
		.get_matches();
    
	// prepare input or quit
	let bam_file    = matches.value_of("BAM").unwrap();
	let bed_file    = matches.value_of("POSITIONS").unwrap();
	let out_file    = matches.value_of("OUT").unwrap();
	let ref_file    = matches.value_of("REF").unwrap_or("NONE");
	let bam_threads = matches.value_of("THREADS").unwrap().parse::<usize>().unwrap();

	if matches.is_present("BAMBAM") {
		bambam::bam_bam_inda_house();
	}
	let positions : HashMap<String,Vec<u64>> = 	parse_bed_file(bed_file);
	let analysis_result = analyze_bam_positions(bam_file,&positions, ref_file, &bam_threads);
	// here we box the error, so in case the writing does
	// not work we get a return of the error from the process back
	if let Err(err) = print_pileup(
			&analysis_result,
			out_file,
			crate_version!(),
			crate_authors!(),&args_string
		){
		eprintln!("{}", err);
        process::exit(1);
    }
	
}
