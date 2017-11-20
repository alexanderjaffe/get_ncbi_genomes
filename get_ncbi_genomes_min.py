'''CHECKS NCBI FOR GENOMES + DOWNLOADS'''

from Bio import Entrez, SeqIO
from bs4 import BeautifulSoup
import time
import lxml
import subprocess
import os
import argparse

def get_genome_ids(search_term):

    # construct + execute entrez search
    handle = Entrez.esearch(db='genome', term=search_term, RetMax=100000)
    result = Entrez.read(handle)
    handle.close()

    return result['IdList']

def get_sequence_ids(genome_id):

    # link genome id -> assembly id
    handle = Entrez.elink(dbfrom = "genome", db="assembly", id=genome_id)
    soup = BeautifulSoup(handle, "xml")
    handle.close()
    # deal with request error
    try: assembly_id = soup.find("Link").find("Id").string
    except: return None, []

    # link assembly to sequence ids
    handle2 = Entrez.elink(dbfrom = "assembly", db="nucleotide", id=assembly_id)
    soup2 = BeautifulSoup(handle2, "xml")
    handle2.close()
    # parse out wgs master and contig names
    for linkset in soup2.find_all("LinkSetDb"):
        set_name = linkset.find("LinkName").string
        # get contig names
        if set_name == "assembly_nuccore_insdc":
            contigs = [link.string for link in linkset.find_all("Id")]
        # get wgs master record
        elif (set_name in ["assembly_nuccore_wgsmaster", "assembly_nuccore_refseq"]):
            wgs_master = linkset.find("Id").string

    return wgs_master, contigs

def get_parse_gbk(sequence_id):

    author_list = []
    # construct/execute gbk request
    handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="gb", retmode="text")
    # get author string
    for record in SeqIO.parse(handle,"genbank"):
        for reference in record.annotations["references"]:
            author_list.append(reference.authors)
        return record, author_list

def download_fasta(sequence_ids, path, basename):

    # construct/execute gbk request
    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta")
    outfile = open(path + "/ncbi_fasta/" + basename + ".fa","w")

    # write out sequences to file
    for seq in SeqIO.parse(handle,"fasta"):
        outfile.write(">" + str(seq.description) + "\n" + str(seq.seq) + "\n")
    handle.close()
    outfile.close()

def record(genome_id, seq_id, gbk, authors, path):

    # take first author of last reference
    author = authors[-1].split(", ")[0]
    source = gbk.annotations["source"]
    #date = time.strftime("%m/%d/%Y")
    # write out to master genome list
    with open(path + "/ncbi_genomes.master","a") as genome_master:
        genome_master.write(genome_id + "\t" + author + "\t" + source + "\n")

'''def process_fasta(path, basename):

    in_path = path + "fasta/" + basename + ".fa"
    # run prodigal
    out_path = path + "prodigal/" + basename
    run_prodigal = "prodigal -q -a " + out_path + ".faa -i " + in_path + " -o " \
        + out_path + ".out"
    subprocess.call(run_prodigal, shell=True)'''

def main():

    # set globs
    __author__ = "Alexander L. Jaffe"
    parser = argparse.ArgumentParser(description='''CHECKS NCBI FOR GENOMES + DOWNLOADS''')
    parser.add_argument('-s','--search_terms', help='Path to search terms file.',required=True)
    parser.add_argument('-o', '--output_directory', help='Path to output directory.', required=True)
    parser.add_argument('-f', '--fasta', help='y/n download genomes.', required=True)
    args = parser.parse_args()
    indir = args.search_terms
    outdir = args.output_directory

    Entrez.email = input("Enter your email address: ")
    downloaded_count = 0
    if args.fasta == "y":
        os.makedirs(outdir + "/ncbi_fasta")
    os.makedirs(outdir + "/ncbi_gbk")

    # define search term list + results
    with open(args.search_terms) as infile:
        lines = infile.readlines()
    search_terms = [line.strip("\n") for line in lines]

    search_results = []
    # get genome ids for each search term
    for term in search_terms:
        search_results += get_genome_ids(term)
        # delay next request
        time.sleep(1)

    # account for overlapping search terms
    unique_results = list(set(search_results))

    # process each genome id
    for genome_id in unique_results:

        print("Trying new genome %s..." %(genome_id))
        # get corresponding sequence ids
        try: wgs_master, contigs = get_sequence_ids(genome_id)
        except:
            print("Network error on %s" %(genome_id))
            continue
        # first check authorship
        if (wgs_master) and (len(contigs)>0):
            # get genbank record + author string
            try: gbk, author_list = get_parse_gbk(wgs_master)
            except:
                print("Network error on %s" %(genome_id))
                continue
            # concatenate author list
            author_string = " ".join(list(set(author_list)))
            # delay next request
            time.sleep(1)

            print("Writing Genbank file for genome %s." %(genome_id))
            # write out genbank file
            basename = "genome_" + genome_id
            outfile = open(outdir + "/ncbi_gbk/" + basename + ".gbk","w")
            SeqIO.write(gbk, outfile, "genbank")

            if args.fasta == "y":
                print("Downloading %d contigs for WGS record %s (genome %s)." \
                    %(len(contigs), str(wgs_master), str(genome_id)))
                # write out fasta
                try: download_fasta(contigs, outdir, basename)
                except:
                    print("Network error on %s" %(genome_id))
                    continue
            # record to master lists
            record(genome_id, wgs_master, gbk, author_list, outdir)
            # predict proteins
            #process_fasta(indir, basename)
            # increment DL count
            downloaded_count+=1

        # delay next request
        time.sleep(1)

    # notify me
    print(str(downloaded_count) + " genomes downloaded from NCBI.")

if __name__ == '__main__':
	main()
