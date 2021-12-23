# AgnosticBamCounter (abc)


The goal was to have an agnostic examination of positions as I found e.g. PacBio HIFI read analysis
of `bamtools mpileup` to be sometimes questionable.

Emanuel Schmid-Siegert    
The Agnostic Bam Counter (abc) determines at a given positions the count of observed nucleotides. 
Simply supply a bed file with positions and obtains a tsv file with counts for each ATCG and reference if provided

```
abc 0.3.1
Emanuel Schmid-Siegert
The Agnostic Bam Counter (abc) determines at a given positions the count of observed nucleotides.
Simply supply a bed file with positions and obtains a tsv file with counts for each ATCG and reference if provided.
If a reference file is provided it evaluates further the reference and variant allele frequency.
This is though always a sum over all potentially multi-allelic site.

Note: 
 - one needs to supply single nucleotide positions and not a range
 - it will sort both, by chromosome and then position - this might pose problems with funky contig names!
 - it is similar to BED 0 based so a SAM position 1000 would translate to 999-1000


USAGE:
    abc [OPTIONS] --bam <BAM> --outfile <TSV> --positions <BED>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -b, --bam <BAM>            aligned short or long reads
    -o, --outfile <TSV>        the output file which will be in tsv format
    -p, --positions <BED>      A bed file with positions
    -r, --reference <FASTA>    if reference in fasta format provided, reference nucleotide 
                               is provided for each positions, otherwise NA
    -t, --threads <INT>        number of threads used in thread-pool for querying positions  [default: 1]

```

## Limitations

 - the bed-file currently accepts only single nucleotide positions, no ranges.
 - systematic testing missing so far

## Output

Example:

```
##	abc:	0.2.2											
##	author:	Emanuel Schmid-Siegert											
##	date:	Sat, 11 Sep 2021 07:35:03 +0200										
##	command:	target/release/abc --bam test/test_mpileup.bam --positions test/test_mpileup.bed --outfile test/test_mpileup.out.tsv --reference test/test_ref.fa									
# chromosome	start	end	reference	A	T	C	G	ambigious	ins	del	depthVAF	RAF
17	301	302	T	1	8	0	0	0	8	0	17	0.5294	0.4706
17	302	303	G	0	0	0	17	0	0	0	17	0.0000	1.0000
17	827	828	T	0	2	11	0	0	0	0	13	0.8462	0.1538
17	1868	1869	A	11	7	0	0	0	0	0	18	0.3889	0.6111
17	2040	2041	G	13	0	0	10	0	0	0	23	0.5652	0.4348
21	10402805	10402806	A	87	0	0	0	0	0	0	87	0.00001.0000
21	10402975	10402976	C	2	2	333	0	0	0	0	337	0.01190.9881

```

The header contains the version and author of the program as well as the execution time and the used command.
We have then in the table the following columns:

 - contig/chromosome name
 - the start (0 based)
 - the end (0 based)
 - the reference nucleotide if reference was provided, otherwise "NA"
 - the count of As
 - the count of Ts
 - the count of Cs
 - the count of Gs
 - the count of ambigious based, e.g. N or X
 - number of reads with an insertion
 - number of reads with a deletion at that position
 - the depth at that position
 - information if mutated or not (true only possible with reference)
