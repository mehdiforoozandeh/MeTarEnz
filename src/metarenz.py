import multiprocessing
multiprocessing.freeze_support()
import warnings, sys, os
warnings.filterwarnings("ignore")
import assembler, screener, blasters, PredNGn
import numpy as np
import pandas as pd


def assembly_to_screening(input_reads_file, screening_database, output_name=None, sra_or_fq='sra', BitScore_filter=50, 
                            Quality_control=True, number_of_threads=1, minimum_contig_length=300, lipase_pred=False):
    
    if output_name ==None:
        if sra_or_fq == 'sra':
            output_name = str(input_reads_file)
            if '/' in output_name:
                output_name = output_name.split('/')[-1]
            elif '//' in output_name:
                output_name = output_name.split('//')[-1]
        elif sra_or_fq == 'fq':
            output_name = str(input_reads_file[0])
            if '/' in output_name:
                output_name = output_name.split('/')[-1]
            elif '//' in output_name:
                output_name = output_name.split('//')[-1]
            
    print('Assembly...')
    if sra_or_fq == 'sra':  
        # in this case input_reads_file should be one SRA file (w/ or without suffix)
        assembler.from_sra_to_cont(input_reads_file, out_name=output_name, QC=Quality_control, 
                                   n_threads=number_of_threads, min_cont_len=minimum_contig_length)
        
    elif sra_or_fq == 'fq': 
        # in this case input_reads_file should be a list of 1 or 2 strings(names or addresses of fastq files)
        assembler.from_fq_to_cont(input_reads_file, out_name=output_name, QC=Quality_control, 
                                  n_threads=number_of_threads, min_cont_len=minimum_contig_length) 
        
    
    assembler.collect_all_res(output_name)
    contig_file_dir = 'MeTarEnz_Results/' + output_name + '/assembly_results/final.contigs.fa'
    
    print('Screening...')
    blasters.makedb(screening_database)
    screening_results = screener.multi_thread_screener(contig_file_dir, BitScore_filter, n_threads=number_of_threads, 
                                                       seq_type='nuc', max_n_lines=2*10**5)
    
    if lipase_pred == True and screening_results.shape[0] > 0:
        lip_pred_results = []
        print('Predicting Lipase Optima...')
        for i in range(screening_results.shape[0]):
            translated_seq = screening_results['translation'][i]
            Topt, pHopt = PredNGn.single_prediction(translated_seq)
            lip_pred_results.append([Topt, pHopt])
        lip_pred_results = pd.DataFrame(np.reshape(np.array(lip_pred_results), (len(lip_pred_results),2)), 
                                        columns=['lipase optimum temperature', 'lipase optimum pH'])
        screening_results = pd.concat([screening_results, lip_pred_results], axis=1)
    
    if screening_results.shape[0] > 0:
        if not os.path.exists('MeTarEnz_Results/' + output_name):
            os.mkdir('MeTarEnz_Results/' + output_name)
        screening_results.to_csv('MeTarEnz_Results/' + output_name + '/screening_results.csv')
        print('Done! Your screening results have been written on MeTarEnz_Results/{}/screening_results.csv'.format(output_name))
    return screening_results
        
        
def contig_screening(input_contigs_file, screening_database, output_name=None, BitScore_filter=50, number_of_threads=1, lipase_pred=False):
    if output_name ==None:
        output_name = str(input_contigs_file)
        if '/' in output_name:
            output_name = output_name.split('/')[-1]
        elif '//' in output_name:
            output_name = output_name.split('//')[-1]
    
    print('Screening...')
    blasters.makedb(screening_database)
    screening_results = screener.multi_thread_screener(input_contigs_file, BitScore_filter, n_threads=number_of_threads, 
                                                       seq_type='nuc', max_n_lines=2*10**5)
    
    if lipase_pred == True and screening_results.shape[0] > 0:
        lip_pred_results = []
        print('Predicting Lipase Optima...')
        for i in range(screening_results.shape[0]):
            translated_seq = screening_results['translation'][i]
            Topt, pHopt = PredNGn.single_prediction(translated_seq)
            lip_pred_results.append([Topt, pHopt])
        lip_pred_results = pd.DataFrame(np.reshape(np.array(lip_pred_results), (len(lip_pred_results),2)), 
                                        columns=['lipase optimum temperature', 'lipase optimum pH'])
        screening_results = pd.concat([screening_results, lip_pred_results], axis=1)
    
    if screening_results.shape[0] > 0:
        if not os.path.exists('MeTarEnz_Results/' + output_name):
            os.mkdir('MeTarEnz_Results/' + output_name)
        screening_results.to_csv('MeTarEnz_Results/' + output_name + '/screening_results.csv')
        print('Done! Your screening results have been written on MeTarEnz_Results/{}/screening_results.csv'.format(output_name))
    return screening_results
    
def protein_screening(input_protein_file, screening_database, output_name=None, BitScore_filter=50, number_of_threads=1, lipase_pred=False):
    if output_name ==None:
        output_name = str(input_protein_file)
        if '/' in output_name:
            output_name = output_name.split('/')[-1]
        elif '//' in output_name:
            output_name = output_name.split('//')[-1]
        
    print('Screening...')
    blasters.makedb(screening_database)
    screening_results = screener.multi_thread_screener(input_protein_file, BitScore_filter, n_threads=number_of_threads, 
                                                       seq_type='pep', max_n_lines=2*10**5)
    
    if lipase_pred == True and screening_results.shape[0] > 0:
        lip_pred_results = []
        print('Predicting Lipase Optima...')
        for i in range(screening_results.shape[0]):
            translated_seq = screening_results['query_seq'][i]
            Topt, pHopt = PredNGn.single_prediction(translated_seq)
            lip_pred_results.append([Topt, pHopt])
        lip_pred_results = pd.DataFrame(np.reshape(np.array(lip_pred_results), (len(lip_pred_results),2)), 
                                        columns=['lipase optimum temperature', 'lipase optimum pH'])
        screening_results = pd.concat([screening_results, lip_pred_results], axis=1)
    
    if screening_results.shape[0] > 0:
        if not os.path.exists('MeTarEnz_Results/' + output_name):
            os.mkdir('MeTarEnz_Results/' + output_name)
        screening_results.to_csv('MeTarEnz_Results/' + output_name + '/screening_results.csv')
        print('Done! Your screening results have been written on MeTarEnz_Results/{}/screening_results.csv'.format(output_name))
    return screening_results





def interactive():
    os.system('clear')
    q1 = input('''** Please choose among the following methods (insert the option's number and press Enter):
               
1)Assembly and Screening from SRA file
2)Assembly and Screening from FASTQ file
3)Screening contigs file (*.fasta)
4)Screening protein file (*.fasta)
5)Prediction of Lipases Temperature and pH optima (*.fasta)\n\n''')
    
    try:
        if str(q1) == '1':
            opt1_questions = {'input':'__Input SRA File:\n', 'db':'__Screening Database File (*.fasta):\n', 
                      'bs':'__Bit-Score Filter (Default: 50):\n', 'output':'__Output File Name (Default: MeTarEnz_Results/"input_file_name"):\n',
                      'qc':'__Report Quality Control Results (y/n) (Default: True): \n', 'n_t':'__Number of Threads (Default: 1): \n', 
                      'mcl':'__Minimum Contig Length (Default: 300): \n', 
                      'lipred': '__Predict Lipase temperature and pH optima ,only if your target enzymes are lipases (y/n) (Default: False)? \n'}
            for k,v in opt1_questions.items():
                if k == 'qc':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt1_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt1_questions[k] = False
                elif k == 'lipred':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt1_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt1_questions[k] = False
                else:   
                    opt1_questions[k] = input(v)
                
            assembly_to_screening(input_reads_file=opt1_questions['input'], screening_database=opt1_questions['db'], 
                        output_name=opt1_questions['output'], sra_or_fq='sra', BitScore_filter=int(opt1_questions['bs']), Quality_control=opt1_questions['qc'], 
                        number_of_threads=int(opt1_questions['n_t']), minimum_contig_length=int(opt1_questions['mcl']), lipase_pred=opt1_questions['lipred'])
                
        elif str(q1) == '2':
            opt2_questions = {'input':'__Input FASTQ Files (if paired-end, separate two files with only a comma):\n', 
                      'db':'__Screening Database File (*.fasta):\n', 'bs':'__Bit-Score Filter (Default: 50):\n', 
                      'output':'__Output File Name (Default: MeTarEnz_Results/"input_file_name"):\n', 'qc':'__Report Quality Control Results (y/n) (Default: Yes): \n', 
                      'n_t':'__Number of Threads (Default: 1): \n', 'mcl':'__Minimum Contig Length (Default: 300): \n', 
                      'lipred':'__Predict Lipase temperature and pH optima ,only if your target enzymes are lipases (y/n) (Default: False)? \n'}
            for k,v in opt2_questions.items():
                
                if k == 'input':
                    temp_inp = input(v)
                    if ',' in temp_inp:
                        opt2_questions[k] = temp_inp.split(',')
                    else: 
                        opt2_questions[k] = [temp_inp]
                        
                elif k == 'qc':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt2_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt2_questions[k] = False
                        
                elif k == 'lipred':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt2_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt2_questions[k] = False
                        
                else:   
                    opt2_questions[k] = input(v)

            assembly_to_screening(input_reads_file = opt2_questions['input'], screening_database= opt2_questions['db'], 
                        output_name=opt2_questions['output'], sra_or_fq='fq', BitScore_filter=int(opt2_questions['bs']), Quality_control=opt2_questions['qc'], 
                        number_of_threads=int(opt2_questions['n_t']), minimum_contig_length=int(opt2_questions['mcl']), lipase_pred=opt2_questions['lipred'])
                    
        elif str(q1) == '3':
            opt3_questions = {'input':'__Input Contig File (*.fasta):\n', 
                      'db':'__Screening Database File (*.fasta):\n', 'bs':'__Bit-Score Filter (Default: 50):\n', 
                      'output':'__Output File Name (Default: MeTarEnz_Results/"input_file_name"):\n', 'n_t':'__Number of Threads (Default: 1): \n', 
                      'lipred':'__Predict Lipase temperature and pH optima ,only if your target enzymes are lipases (y/n) (Default: False)? \n'}
            for k,v in opt3_questions.items():
                if k == 'lipred':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt3_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt3_questions[k] = False
                else:
                    opt3_questions[k] = input(v)
            contig_screening(input_contigs_file=opt3_questions['input'], screening_database=opt3_questions['db'], 
                        output_name=opt3_questions['output'], BitScore_filter=int(opt3_questions['bs']), 
                        number_of_threads=int(opt3_questions['n_t']), lipase_pred=opt3_questions['lipred'])
            
        elif str(q1) == '4':
            opt4_questions = {'input':'__Input Protein File (*.fasta):\n', 
                      'db':'__Screening Database File (*.fasta):\n', 'bs':'__Bit-Score Filter (Default: 50):\n', 
                      'output':'__Output File Name (Default: MeTarEnz_Results/"input_file_name"):\n', 'n_t':'__Number of Threads (Default: 1): \n', 
                      'lipred':'__Predict Lipase temperature and pH optima ,only if your target enzymes are lipases (y/n) (Default: False)? \n'}
            for k,v in opt4_questions.items():
                if k == 'lipred':
                    temp_inp = input(v)
                    if temp_inp.lower() == 'y':
                        opt4_questions[k] = True
                    elif temp_inp.lower() == 'n':
                        opt4_questions[k] = False
                else:
                    opt4_questions[k] = input(v)
            protein_screening(input_protein_file=opt4_questions['input'], screening_database=opt4_questions['db'], 
                        output_name=opt4_questions['output'], BitScore_filter=int(opt4_questions['bs']), 
                        number_of_threads=int(opt4_questions['n_t']), lipase_pred=opt4_questions['lipred'])
                
        elif str(q1) == '5':
            opt5_questions = {'input':'__Input Protein File (*.fasta):\n', 
                              'output':'__Output File Name (Default: MeTarEnz_Results/"input_file_name"):\n'}
            for k,v in opt5_questions.items():
                opt5_questions[k] = input(v)
            results = PredNGn.fasta_prediction(opt5_questions['input'])
            results = pd.DataFrame(results)
            results.to_csv('MeTarEnz_Results/' + opt5_questions['output'] + '/lipase_prediction.csv')
    except:
        print('''Error: 
Something went wrong! 
Please check your input files and/or command. 
You can use -h or --help for more information about the usage or read the manual for detailed instructions.
''')

  
if __name__ == '__main__':
    with open('help', 'r') as hm:
        helpme = hm.read()
    
    if len(sys.argv)==1: print('MeTarEnz: missing operand \nTry -h or --help for more information.')
    
    elif len(sys.argv)==2:
        if sys.argv[1]=='-h' or sys.argv[1]=='--help': print(helpme)
        elif sys.argv[1]=='-int' or sys.argv[1]=='--interactive': interactive()
        else:
            print('MeTarEnz: missing operand \nTry -h or --help for more information.')
        
    elif len(sys.argv)>2:  
        try:    
            if sys.argv[1]== 'sras' or sys.argv[1]=='sra_screening':
                
                oblig_args = {'input': None, 'db': None}
                opt_args = {'bs':50, 'output':None, 'lipred':False, 'qc': True, 'n_t':1, 'mcl':300}
                
                if '-i' in sys.argv:
                    oblig_args['input'] = str(sys.argv[sys.argv.index('-i')+1])
                    
                if '-db' in sys.argv:
                    oblig_args['db']= str(sys.argv[sys.argv.index('-db')+1])
                    
                if '-bs' in sys.argv:
                    opt_args['bs']= int(sys.argv[sys.argv.index('-bs')+1])
                    
                if '-lip' in sys.argv:
                    opt_args['lipred']= True
                
                if '-o' in sys.argv:
                    opt_args['output']= str(sys.argv[sys.argv.index('-o')+1])
                
                if '-noqc' in sys.argv:
                    opt_args['qc']= False         
                        
                if '-t' in sys.argv:
                    opt_args['n_t']= int(sys.argv[sys.argv.index('-t')+1])
                    
                if '--min-contig-len' in sys.argv:
                    opt_args['mcl']= int(sys.argv[sys.argv.index('--min-contig-len')+1])
                if '-mcl' in sys.argv:
                    opt_args['mcl']= int(sys.argv[sys.argv.index('-mcl')+1])
                    
                if None in oblig_args.values():
                    print("Missing necessary arguments (input -i or the screening database -db).")
                else:
                    assembly_to_screening(input_reads_file = oblig_args['input'], screening_database= oblig_args['db'], 
                                        output_name=opt_args['output'], sra_or_fq='sra', BitScore_filter=opt_args['bs'], Quality_control=opt_args['qc'], 
                                        number_of_threads=opt_args['n_t'], minimum_contig_length=opt_args['mcl'], lipase_pred=opt_args['lipred'])

            elif sys.argv[1]== 'fqs' or sys.argv[1]=='fastq_screening':
                
                oblig_args = {'input': None, 'db': None, 'end':None}
                opt_args = {'bs':50, 'output':None, 'lipred':False, 'qc': True, 'n_t':1, 'mcl':300}
                
                if '-s' in sys.argv:
                    oblig_args['input'] = [str(sys.argv[sys.argv.index('-s')+1])]
                    oblig_args['end'] = 'single'
                elif '-p' in sys.argv:
                    oblig_args['input']= [str(sys.argv[sys.argv.index('-p')+1]), str(sys.argv[sys.argv.index('-p')+2])]
                    oblig_args['end'] = 'pair'
                
                if '-db' in sys.argv:
                    oblig_args['db']= str(sys.argv[sys.argv.index('-db')+1])
                
                if '-bs' in sys.argv:
                    opt_args['bs']= int(sys.argv[sys.argv.index('-bs')+1])
                    
                if '-lip' in sys.argv:
                    opt_args['lipred']= True
                
                if '-o' in sys.argv:
                    opt_args['output']= str(sys.argv[sys.argv.index('-o')+1])
                
                if '-noqc' in sys.argv:
                    opt_args['qc']= False         
                        
                if '-t' in sys.argv:
                    opt_args['n_t']= int(sys.argv[sys.argv.index('-t')+1])
                    
                if '--min-contig-len' in sys.argv:
                    opt_args['mcl']= int(sys.argv[sys.argv.index('--min-contig-len')+1])
                if '-mcl' in sys.argv:
                    opt_args['mcl']= int(sys.argv[sys.argv.index('-mcl')+1])
                    
                if None in oblig_args.values():
                    print("Missing necessary arguments (single -s or pair -p ended fastq with two input files, screening database -db).") 
                else:
                    assembly_to_screening(input_reads_file = oblig_args['input'], screening_database= oblig_args['db'], 
                                        output_name=opt_args['output'], sra_or_fq='fq', BitScore_filter=opt_args['bs'], Quality_control=opt_args['qc'], 
                                        number_of_threads=opt_args['n_t'], minimum_contig_length=opt_args['mcl'], lipase_pred=opt_args['lipred'])
                
                
            elif sys.argv[1]== 'cs' or sys.argv[1]== 'contig_screening':
                oblig_args = {'input': None, 'db': None}
                opt_args = {'bs':50, 'output':None, 'lipred':False, 'n_t':1}
                
                if '-i' in sys.argv:
                    oblig_args['input'] = str(sys.argv[sys.argv.index('-i')+1])
                
                if '-db' in sys.argv:
                    oblig_args['db']= str(sys.argv[sys.argv.index('-db')+1])
                
                if '-bs' in sys.argv:
                    opt_args['bs']= int(sys.argv[sys.argv.index('-bs')+1])
                    
                if '-lip' in sys.argv:
                    opt_args['lipred']= True
                
                if '-o' in sys.argv:
                    opt_args['output']= str(sys.argv[sys.argv.index('-o')+1])
                        
                if '-t' in sys.argv:
                    opt_args['n_t']= int(sys.argv[sys.argv.index('-t')+1])
                    
                if None in oblig_args.values():
                    print("Missing necessary arguments (input -i or the screening database -db).")
                else:
                    contig_screening(input_contigs_file=oblig_args['input'], screening_database=oblig_args['db'], 
                                        output_name=opt_args['output'], BitScore_filter=opt_args['bs'], 
                                        number_of_threads=opt_args['n_t'], lipase_pred=opt_args['lipred'])
                    
            
            elif sys.argv[1]== 'ps' or sys.argv[1]== 'protein_screening':
                oblig_args = {'input': None, 'db': None}
                opt_args = {'bs':50, 'output':None, 'lipred':False, 'n_t':1}
                
                if '-i' in sys.argv:
                    oblig_args['input'] = str(sys.argv[sys.argv.index('-i')+1])
                
                if '-db' in sys.argv:
                    oblig_args['db']= str(sys.argv[sys.argv.index('-db')+1])
                
                if '-bs' in sys.argv:
                    opt_args['bs']= int(sys.argv[sys.argv.index('-bs')+1])
                    
                if '-lip' in sys.argv:
                    opt_args['lipred']= True
                
                if '-o' in sys.argv:
                    opt_args['output']= str(sys.argv[sys.argv.index('-o')+1])
                        
                if '-t' in sys.argv:
                    opt_args['n_t']= int(sys.argv[sys.argv.index('-t')+1])
                    
                if None in oblig_args.values():
                    print("Missing necessary arguments (input -i or the screening database -db).")
                else:
                    protein_screening(input_protein_file=oblig_args['input'], screening_database=oblig_args['db'], 
                                        output_name=opt_args['output'], BitScore_filter=opt_args['bs'], 
                                        number_of_threads=opt_args['n_t'], lipase_pred=opt_args['lipred'])
            
            elif sys.argv[1]== 'lp' or sys.argv[1]== 'lipase_prediction':
                oblig_args = {'input': None, 'output': None}
                
                if '-i' in sys.argv:
                        oblig_args['input'] = str(sys.argv[sys.argv.index('-i')+1])
                if '-o' in sys.argv:
                    oblig_args['output']= str(sys.argv[sys.argv.index('-o')+1])

                if None in oblig_args.values():
                    print("Missing necessary arguments (input -i or the output file name -o).")
                else:
                    results = PredNGn.fasta_prediction(oblig_args['input'])
                    results = pd.DataFrame(results)
                    if not os.path.exists('MeTarEnz_Results/'+oblig_args['output']):
                        os.mkdir('MeTarEnz_Results/'+oblig_args['output'])
                    results.to_csv('MeTarEnz_Results/' + oblig_args['output'] + '/lipase_prediction_results.csv')
                    print('Done! Your lipase prediction results have been written on MeTarEnz_Results/{}/lipase_prediction_results.csv'.format(oblig_args['output']))

            else:
                print('Invalid command!')  
        except:
            print('''Error: 
Something went wrong! Please check your input files and/or command.
You can use -h or --help for more information about the usage or read the manual for detailed instructions. 
Interactive execution of MeTarEnz is possible with -int or --interactive options.''')
    