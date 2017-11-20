# get_ncbi_genomes_min.py
Download genomes from NCBI Genome database using compound queries.

```
usage: get_ncbi_genomes_min.py [-h] -s SEARCH_TERMS -o OUTPUT_DIRECTORY -f FASTA

optional arguments:
  -h, --help            show this help message and exit
  -s SEARCH_TERMS, --search_terms SEARCH_TERMS
                        Path to search terms file.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output directory.
  -f FASTA, --fasta FASTA
                        y/n download genomes.
```

<b>Output:</b>

`ncbi_gbk` : directory containing genbank files for downloaded genomes <br>
`ncbi_fasta` : directory containing fasta files for downloaded genomes <br>
`ncbi_genomes.master` : log file of format genome_id\tauthor\tgenome_name <br>

<b>example</br> `search_terms.txt`:

```
Klebsiella NOT phage
Chloroflexi AND 2015
Aenigmarchaeota
...
```
<b>NB: Recommended to manually check results of search terms before running script. </b>
