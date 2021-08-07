# AgnosticBamCounter (abc)

**this tool is very beta** !


The goal was to have an agnostic examination of positions as I found e.g. PacBio HIFI read analysis
of `bamtools mpileup` to be sometimes questionable.

abc 0.1.0    
Emanuel Schmid-Siegert    
The Agnostic Bam Counter (abc) determines at a given positions the count of observed nucleotides. 
Simply supply a bed file with positions and obtains a tsv file with counts for each ATCG and reference if provided

```
USAGE:
    abc [FLAGS] --bam <BAM> --outfile <TSV> --positions <BED>

FLAGS:
    -r, --reference    if reference in fasta format provided, reference nucleotide 
                       is provided for each positions, otherwise NA
    -h, --help         Prints help information
    -V, --version      Prints version information

OPTIONS:
    -b, --bam <BAM>          aligned short or long reads
    -o, --outfile <TSV>      the output file which will be in tsv format
    -p, --positions <BED>    A bed file with positions
```

## Limitations

 - the bed-file currently accepts only single nucleotide positions, no ranges.
 - not tested for indels
 - systematic testing missing so far

