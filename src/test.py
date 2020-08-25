print ('Running the tests to see if MeTarEzn Actually Works ...')
import os
from colored import fg, bg, attr
import shutil  

def inpy_test():
    try:
        import metarenz, PredNGn
        import pandas as pd
    except:
        print('Error:   Failed at import level')
        exit()

    c = fg('green')
    r = attr('reset')
    try:
        print (c + 'testing SRA assembly_to_screening...' + r)
        metarenz.assembly_to_screening(input_reads_file = 'test_inputs/SRR12198400', screening_database= 'test_inputs/testseqs/test_db.fasta', 
                        output_name='test', sra_or_fq='sra', BitScore_filter=50, 
                        Quality_control=True, number_of_threads=1, lipase_pred=True)

        print (c + 'testing FQ assembly_to_screening...' + r)
        metarenz.assembly_to_screening(input_reads_file = ['test_inputs/p_file1_1.fq', 'test_inputs/p_file1_2.fq'], 
                        screening_database= 'test_inputs/testseqs/test_db.fasta', 
                        output_name='test', sra_or_fq='fq', BitScore_filter=50, 
                        Quality_control=True, number_of_threads=1, lipase_pred=True)

        print (c + 'testing contig_screening...' + r)
        metarenz.contig_screening(input_contigs_file='test_inputs/testseqs/test_nuc.fasta', screening_database='test_inputs/testseqs/test_db.fasta', 
                        output_name='test', BitScore_filter=700, 
                        number_of_threads=1, lipase_pred=True)

        print (c + 'testing protein_screening...' + r)
        metarenz.protein_screening(input_protein_file='test_inputs/testseqs/test_aa.fasta', screening_database='test_inputs/testseqs/test_db.fasta', 
                        output_name='test', BitScore_filter=800, 
                        number_of_threads=1, lipase_pred=True)

        print (c + 'testing lipase prediction...' + r)
        results = PredNGn.fasta_prediction('test_inputs/testseqs/test_lip.fasta')
        results = pd.DataFrame(results)
        if not os.path.exists('MeTarEnz_Results/test'):
            os.mkdir('MeTarEnz_Results/test')
        results.to_csv('MeTarEnz_Results/test/lipase_prediction_results.csv')

        if os.path.exists('MeTarEnz_Results/test'):
            shutil.rmtree('MeTarEnz_Results/test')
        print(c + 'Done! \nWorks like a charm!' + r)

    except:
        if os.path.exists('MeTarEnz_Results/test'):
            shutil.rmtree('MeTarEnz_Results/test')
        print(fg('red') + 'The test did not work out well!' + r)
        exit()

def inshell_test():
    c = fg('green')
    r = attr('reset')
    try:
        print (c + 'testing SRA assembly_to_screening...' + r)
        os.system('python3 metarenz.py sras -i test_inputs/SRR12198400 -db test_inputs/testseqs/test_db.fasta -o test -bs 100 -lip')

        print (c + 'testing FQ assembly_to_screening...' + r)
        os.system('python3 metarenz.py fqs -p test_inputs/p_file1_1.fq test_inputs/p_file1_2.fq -db test_inputs/testseqs/test_db.fasta -o test -bs 100 -lip')

        print (c + 'testing contig_screening...' + r)
        os.system('python3 metarenz.py cs -i test_inputs/testseqs/test_nuc.fasta -db test_inputs/testseqs/test_db.fasta -o test -bs 700 -lip')

        print (c + 'testing protein_screening...' + r)
        os.system('python3 metarenz.py ps -i test_inputs/testseqs/test_aa.fasta -db test_inputs/testseqs/test_db.fasta -o test -bs 800 -lip')

        print (c + 'testing lipase prediction...' + r)
        os.system('python3 metarenz.py lp -i test_inputs/testseqs/test_lip.fasta -o test')

        if os.path.exists('MeTarEnz_Results/test'):
            shutil.rmtree('MeTarEnz_Results/test')
        print(c + 'Done! \nWorks like a charm!' + r)

    except:
        if os.path.exists('MeTarEnz_Results/test'):
            shutil.rmtree('MeTarEnz_Results/test')
        print(fg('red') + 'The test did not work out well!' + r)
        exit()

if __name__ == '__main__':
        inshell_test()
