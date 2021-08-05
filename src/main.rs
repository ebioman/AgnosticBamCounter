use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
//use bio::io::{bed::Reader};
use std::collections::HashMap;
use std::path::Path;
use rust_htslib::{bam, bam::Read};
use std::str::from_utf8;



#[derive(Debug,Default)]
/// structure which holds our
/// counts of nucleotides for 
/// any given genomic position
struct Pileup {
    /// chromosome
    chromosome: String,
    /// position
    position: u64,
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
	/// depth at that given position
	depth: u32,
}

// takes the generates pileups in our specific
// format and prints them according to our wishes
// currently only tsv, might add more possibilities
// in the future
fn print_pileup(result:Vec<Pileup>){
	if result.len() == 0 {
		eprintln!("INFO: no pileup information aquired - good bye")
	}
	for element in result {
		println!{"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", element.chromosome,element.position,element.reference,element.nuc_a,element.nuc_t,element.nuc_c,element.nuc_g,element.ambigious,element.depth}
	}
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
		if &record.end()-&record.start()>1 {
			panic!("ERROR: entry {:?} contained not a position but a range!",&record);
		}
		bed_result
			.entry(record.chrom().to_string())
			.or_insert_with(Vec::new)
			.push(record.start());
	
	};
	bed_result
}

// This reads now the BAM file and uses the positions
// to check for each position the observed number of nucleotides
fn analyze_bam_positions (input: &str , positions: &HashMap<String,Vec<u64>>) -> Vec<Pileup>{
	assert!(
		Path::new(input).exists(),
		"ERROR: BAM file {} does not exist!",
		input,
	);
	let mut overview : Vec<Pileup> = Vec::new();
	let mut bam_file = bam::IndexedReader::from_path(input).unwrap();
	for (chr,collection) in positions {
		for pos in collection {
			let tmp_result = fetch_position(&mut bam_file,chr,pos);
			overview.push(tmp_result);
		}
	}
	overview
}

// this sub-function gets the iterator result from the positions
// and return the result of pileup analysis at that given position
// allows to potentially already parallelize over each position analyzed
fn fetch_position(bam: &mut bam::IndexedReader, chr: &str, pos:&u64) -> Pileup {
	// NOTE:
	// SAM is 1 based and not 0 based, need to correct for that in
	let start = *pos -1 ;
	let end   = *pos ;
	// this obtains now the pileup at that
	// given position
	bam.fetch((chr,start,end)).expect("ERROR: could not fetch region");
	// currently we ignore completely clipping
	let mut collection : Pileup = Default::default();
	collection.chromosome = chr.to_string();
	collection.position   = *pos;
	for pile in bam.pileup().map(|x| x.expect("ERROR: could not parse BAM file")){
		// now we only care about the position we inquire
		if pile.pos() as u64 == start {
			collection.depth = pile.depth();
			for alignment in pile.alignments(){
				// sometimes we get reads with 0 length, no idea why
				if alignment.record().seq_len() == 0 { continue }
				// some have none, no idea why
				
				let qpos = match  alignment.qpos(){
					Some(q) => q,
					_ =>  continue,
				};
				let nuc  = &alignment.record().seq().as_bytes()[qpos..qpos+1];
				let nuc1 = from_utf8(nuc).unwrap();
				match nuc1 {
					"A" => collection.nuc_a = collection.nuc_a +1  ,
					"T" => collection.nuc_t = collection.nuc_t +1  ,
					"C" => collection.nuc_c = collection.nuc_c +1  ,
					"G" => collection.nuc_g = collection.nuc_g +1  ,
					_ => collection.ambigious = collection.ambigious +1  ,
				}
			}
		}
	};
	collection

}

fn main() {
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
			.takes_value(false)
            .required(false))
		.get_matches();
    
	// prepare input or quit
	let bam_input = matches.value_of("BAM").unwrap();
	let bed_file  = matches.value_of("POSITIONS").unwrap();
	let _out_file  = matches.value_of("OUT").unwrap();
	let _ref_file  = matches.value_of("REF").unwrap_or("NONE");
	
	let positions : HashMap<String,Vec<u64>> = 	parse_bed_file(&bed_file);
	let analysis_result = analyze_bam_positions(&bam_input,&positions);
	print_pileup(analysis_result);
}
