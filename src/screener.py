import blasters, funcs
import multiprocessing as mp
import numpy as np
import pandas as pd
from functools import partial
import random, string, os, shutil

import Bio.SeqIO as SeqIO

def screen_nuc(filename, bs_filter):

    blastresname = 'temp/result_blastx_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7)) + '.txt'
    blasters.blastx(filename, 'temp/db.fasta', blastresname)
    df = blasters.read_blast_results(blastresname)
    
    matches = []
    input_seqs = funcs.read_fasta(filename)
    sixpack = funcs.gen_frames(input_seqs)
    
    # pruning non-homologous results
    # ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/

    to_drop = []
    for i in range(len(df)):
        if float(df.evalue[i]) > 1e-6:
            to_drop.append(i)
        elif float(df.bitscore[i]) < 50:
            to_drop.append(i)
    df = df.drop(to_drop)
    df = df.reset_index()

    # choose the best match
    tlis = {}
    for i in range(len(df)):
        if str(df.qseqid[i]) not in tlis.keys():
            tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]
        else:
            if tlis[str(df.qseqid[i])][1] < float(df.bitscore[i]):
                tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]
    p = list(np.arange(len(df)))
    q = [bs[0] for bs in tlis.values()]
    drop_list = [item for item in p if item not in q]
    df = df.drop(drop_list)
    df = df.reset_index()
    df = df.drop(['level_0', 'index'], axis=1)
    
    # Filter Based on Bit-Score
    for i in range(len(df)):
        qid = df.qseqid[i]
        qseq = df.qseq[i].replace('-', '')
        for frame, trans in sixpack[qid].items():
            trans = trans.split('*')
            for contig in trans:
                if qseq[:20] in contig:
                    match_start = contig.find(qseq[:20])
                    bef = contig[:match_start]
                    if 'M' in bef:
                        bef = bef[::-1]
                        methionin = bef.find('M')
                        translated_region = contig[(match_start - methionin - 1):]
                    else:
                        translated_region = contig[contig.find(qseq[:20]):]
                    if float(df.bitscore[i]) > bs_filter:
                        matches.append([qid, translated_region, frame, df.sseqid[i], df.bitscore[i], df.evalue[i]])
    matches = pd.DataFrame(matches, columns=['query_id', 'translation', 'frame', 'source_seq_id' ,'bitscore', 'evalue'])
    return matches

def screen_pep(filename, bs_filter):
    
    blastresname = 'temp/result_blastp_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7)) + '.txt'
    blasters.blastp(filename, 'temp/db.fasta', blastresname)
    df = blasters.read_blast_results(blastresname)
    
    # pruning non-homologous results
    # ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/

    matches = []
    
    to_drop = []
    for i in range(len(df)):
        if float(df.evalue[i]) > 1e-6:
            to_drop.append(i)
        elif float(df.bitscore[i]) < 50:
            to_drop.append(i)
    df = df.drop(to_drop)
    df = df.reset_index()
    
    # choose the best match
    tlis = {}
    for i in range(len(df)):
        if str(df.qseqid[i]) not in tlis.keys():
            tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]
        else:
            if tlis[str(df.qseqid[i])][1] < float(df.bitscore[i]):
                tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]
    p = list(np.arange(len(df)))
    q = [bs[0] for bs in tlis.values()]
    drop_list = [item for item in p if item not in q]
    df = df.drop(drop_list)
    df = df.reset_index()
    df = df.drop(['level_0', 'index'], axis=1)
    
    # Filter Based on Bit-Score
    for i in range(len(df)):
        qid = df.qseqid[i]
        qseq = df.qseq[i].replace('-', '')
        
        if float(df.bitscore[i]) > bs_filter:
            matches.append([qid, qseq, df.sseqid[i], df.bitscore[i], df.evalue[i]])
        
    if len(matches)!= 0:
        matches = pd.DataFrame(matches, columns=['query_id', 'query_seq', 'source_seq_id' ,'bitscore', 'evalue'])
    return matches


def multi_thread_screener(filename, bs_filter, n_threads=1, seq_type='nuc', max_n_lines=2*10**5):
    tempor_fol_name = 'temp/temp_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7))
    try:
        
        os.mkdir(tempor_fol_name)
        with open(filename) as fin:
            fout = open(tempor_fol_name + '/tempgene0.fasta', "w")
            len_fout = 0
            short_memory = []
            for i, line in enumerate(fin):
                
                if line[0] == '>' and len(short_memory) == 0:
                    short_memory.append(line)
                elif line[0] != '>' and len(short_memory) > 0:
                    short_memory.append(line)
                elif line[0] == '>' and len(short_memory) > 0:
                    for ll in short_memory:
                        fout.write(ll)
                        len_fout +=1
                    short_memory = []
                    short_memory.append(line)
                    
                if len_fout >= max_n_lines:
                    fout.close()
                    fout = open(tempor_fol_name + '/tempgene%d.fasta' % (i / max_n_lines + 1), "w")     
                    len_fout = 0
                    
                    
            fout.close()
        
        temp_files = os.listdir(tempor_fol_name)
        for i in range(len(temp_files)):
            temp_files[i] = tempor_fol_name + '/' + str(temp_files[i])

        if n_threads != 1 and isinstance(n_threads, type(4)) and n_threads > 1:
            if n_threads > mp.cpu_count():
                n_threads = mp.cpu_count()

        pool = mp.Pool(n_threads)
        results = []
        if seq_type == 'nuc':
            N = pool.map(partial(screen_nuc, bs_filter=bs_filter), temp_files)
        elif seq_type == 'pep':
            N = pool.map(partial(screen_pep, bs_filter=bs_filter), temp_files)    
        for i in N:
            results.append(i) 
        
        shutil.rmtree(tempor_fol_name)
        results = pd.concat(results, axis=0)
        results = results.reset_index()
        results = results.drop(['index'], axis=1)
        print('Screening Results Dataframe Shape: '+ str(results.shape))
        return results
    
    except:
        print("Something Went Wrong!")
        shutil.rmtree(tempor_fol_name)
        exit()

    
if __name__ == '__main__':
    from datetime import datetime
    t0 = datetime.now()


    blasters.makedb('testseqs/final_ref.fasta')
    
    for t in range(1,5):
        t0 = datetime.now()
        res1 = multi_thread_screener('testseqs/EGL.fasta', 300, n_threads=t, seq_type='nuc')
        print('Run Time with {} threads = '.format(str(t)), datetime.now() - t0, res1.shape)
    
    for t in range(1,5):
        t0 = datetime.now()
        res1 = multi_thread_screener('testseqs/3214_prot.fasta', 300, n_threads=t, seq_type='pep')
        print('Run Time with {} threads = '.format(str(t)), datetime.now() - t0, res1.shape)

    