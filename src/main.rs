use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::collections::HashMap;
use std::process;
use std::env;
extern crate bambam;

use crate::libs::{print_pileup,parse_bed_file,analyze_bam_positions };
mod libs;

fn main() {
	env_logger::init();

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
