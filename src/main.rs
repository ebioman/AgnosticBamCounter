use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
//use bio::io::{bed::Reader};
use std::collections::HashMap;
use std::path::Path;
use rust_htslib::{bam, bam::Read};
use std::str::from_utf8;
use bio::io::fasta::IndexedReader;
use chrono::{DateTime, Local};
use csv;
use std::error;
use std::process;
use std::env;
use itertools::Itertools;


// Change the alias to `Box<error::Error>`.
type Result<T> = std::result::Result<T, Box<dyn error::Error>>;




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
fn print_pileup(result:&Vec<Pileup>,out: &str, version: &str, author: &str , command: &str )-> Result<i32> {
	let now: DateTime<Local> = Local::now();
	let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(&out)?;
	writer.write_record(&["##","abc:",version,"","","","","",""])?;
    writer.write_record(&["##","author:",author,"","","","","",""])?;
    writer.write_record(&["##","date:",&now.to_rfc2822(),"","","","","",""])?;
    writer.write_record(&["##","command:", command,"","","","","",""])?;    
    writer.write_record(&["# chromosome","position","reference","A","T","C","G","ambigious","depth"])?;
    for element in result.iter() {
        writer.write_record(&[
			element.chromosome.clone(),
			element.position.to_string(),
			element.reference.clone(),
			element.nuc_a.to_string(),
			element.nuc_t.to_string(),
			element.nuc_c.to_string(),
			element.nuc_g.to_string(),
			element.ambigious.to_string(),
			element.depth.to_string()
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
// It sorts as well again by chromosome and position to be safe
fn analyze_bam_positions (input: &str , positions: &HashMap<String,Vec<u64>>, ref_fasta: &str, threads: &u32) -> Vec<Pileup>{
	assert!(
		Path::new(input).exists(),
		"ERROR: BAM file {} does not exist!",
		input,
	);
	let mut overview : Vec<Pileup> = Vec::new();
	let mut bam_file = bam::IndexedReader::from_path(input).unwrap();
	// generate a thread pool for the multithreaded reading and writing of BAM files
	eprintln!("INFO: Parsing BAM file with {} threads",threads);
    let pool = rust_htslib::tpool::ThreadPool::new(*threads).unwrap();
    bam_file.set_thread_pool(&pool).unwrap();
	// this is really not great, I want to open if available the
	// index file and provide to the fetch postion but I cant manage
	// to specify the type in the function definition
	// therefore I open close very often.....
	for chr in positions.keys().sorted() {
		for pos in positions.get(chr).unwrap().iter().sorted(){
			let tmp_result = match ref_fasta { 
			"NONE" => fetch_position(&mut bam_file,chr,pos,None),
			 _     => fetch_position(&mut bam_file,chr,pos,Some(ref_fasta)),
			};
			overview.push(tmp_result);
		}
	}
	overview
}

// this sub-function gets the iterator result from the positions
// and return the result of pileup analysis at that given position
// allows to potentially already parallelize over each position analyzed
fn fetch_position(bam: &mut bam::IndexedReader, chr: &str, pos:&u64, ref_file: Option<&str>) -> Pileup {
	// NOTE:
	// SAM is 1 based and not 0 based, need to correct for that in
	let start = *pos -1 ;
	let end   = *pos ;
	// this obtains now the pileup at that
	// given position
	bam.fetch((chr,start,end)).expect("ERROR: could not fetch region");
	// currently we ignore completely clipping
	let mut collection : Pileup = Default::default();
	if ref_file.is_some() {
		let mut faidx = IndexedReader::from_file(&ref_file.unwrap()).unwrap();
		faidx.fetch(chr,start,end).expect("ERROR: could not fetch interval on reference ");
		let mut ref_seq = Vec::new();
		faidx.read(&mut ref_seq).expect("ERROR: could not read sequence from reference");
		collection.reference = String::from(from_utf8(&ref_seq).unwrap());
	}else{
		collection.reference = String::from("NA");
	}
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
			.takes_value(true)
			.required(false))
		.get_matches();
    
	// prepare input or quit
	let bam_file = matches.value_of("BAM").unwrap();
	let bed_file  = matches.value_of("POSITIONS").unwrap();
	let out_file  = matches.value_of("OUT").unwrap();
	let ref_file  = matches.value_of("REF").unwrap_or("NONE");
	let bam_threads   = matches.value_of("THREADS").unwrap_or("1").parse::<u32>().unwrap();
	let positions : HashMap<String,Vec<u64>> = 	parse_bed_file(&bed_file);
	
	let analysis_result = analyze_bam_positions(&bam_file,&positions, &ref_file, &bam_threads);
	// here we box the error, so in case the writing does
	// not work we get a return of the error from the process back
	if let Err(err) = print_pileup(&analysis_result,&out_file,crate_version!(),crate_authors!(),&args_string){
		eprintln!("{}", err);
        process::exit(1);
    }
	
}
