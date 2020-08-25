import os
import pandas as pd


info = """
Info:

This module contains different functions to interact 
with NCBI's BLASTp and BLASTx tools + BLAST DB maker.
"""


def makedb(fastaname):
	os.system('./makeblastdb -dbtype prot -in '+ str(fastaname)+ ' -out temp/db.fasta') 
 
def blastx(queryfile, database, outfile):
    blast_options = './blastx -outfmt "10 qseqid qlen length sseqid qseq pident evalue bitscore"'
    command = str(blast_options + ' -query ' + queryfile + ' -db ' + database + ' -out ' + outfile)
    os.system(command)
    
def blastp(queryfile, database, outfile):
    blast_options = './blastp -outfmt "10 qseqid qlen length sseqid qseq pident evalue bitscore"'
    command = str(blast_options + ' -query ' + queryfile + ' -db ' + database + ' -out ' + outfile)
    os.system(command)
    
def read_blast_results(blastresname):
    with open(blastresname, 'r') as file:
        blast_res = file.readlines()
        br = []
        for i in blast_res:
            br.append(i[:-1].split(','))
    os.remove(blastresname)
    columns = ['qseqid', 'qlen', 'length', 'sseqid', 
               'qseq', 'pident', 'evalue', 'bitscore']
    df = pd.DataFrame(br, columns=columns)
    return df


if __name__ == '__main__':
    print(info)
    # makedb('testseqs/final_ref.fasta')
    # blastp('testseqs/3214_prot.fasta', 'db.fasta', '3214_prot_blastres.txt')
    # print(read_blast_results('3214_prot_blastres.txt'))