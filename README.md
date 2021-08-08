# AgnosticBamCounter (abc)

**this tool is very beta** !


The goal was to have an agnostic examination of positions as I found e.g. PacBio HIFI read analysis
of `bamtools mpileup` to be sometimes questionable.

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

## Output

Example:

```
##      abc:    0.2.0                                           
##      author: Emanuel Schmid-Siegert                                          
##      date:   Sun, 08 Aug 2021 08:47:45 +0200                                         
##      command:        target/debug/abc --bam test/test_mpileup.bam --outfile test/test_mpileup.out.tsv --positions test/test_mpileup.bed --reference test/test_ref.fa                                         
# chromosome    position        reference       A       T       C       G       ambigious       depth
17      302     T       1       16      0       0       0       17
17      303     G       0       0       0       17      0       17
17      1869    A       11      7       0       0       0       18
17      1870    C       0       0       18      0       0       18
17      1884    G       0       0       0       21      0       21
21      10402985        G       3       0       0       341     0       344
21      10405200        T       0       29      30      0       0       59
```

The header contains the version and author of the program as well as the execution time and the used command.
We have then in the table the following columns:

 - contig/chromosome name
 - the position 0 based
 - the reference nucleotide if reference was provided, otherwise "NA"
 - the count of As
 - the count of Ts
 - the count of Cs
 - the count of Gs
 - the count of ambigious based, e.g. N or X
 - the depth at that position
