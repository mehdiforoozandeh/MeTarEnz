import sys
import os
import shutil
import subprocess
cdp=os.getcwd()

def sra2fq(inp, outp):
    if not os.path.exists('temp/FQs'):
        os.mkdir('temp/FQs')
    
    if not os.path.exists('temp/FQs/'+outp):
        os.mkdir('temp/FQs/'+ outp)
        
    outp = 'temp/FQs/' + outp
    command = "sratk/fastq-dump -I --split-files -O %s %s" % (outp, inp)
    subprocess.check_call(command, shell=True)

def QC_report(inp, outp):
    if not os.path.exists('temp/QC_reports'):
        os.mkdir('temp/QC_reports')
    
    if not os.path.exists('temp/QC_reports/'+outp):
        os.mkdir('temp/QC_reports/'+ outp)
        
    outp = 'temp/QC_reports/' + outp
    command = "FastQC/fastqc -o %s -f fastq %s" % (outp, inp)
    subprocess.check_call(command, shell=True)
    
def single_assembly(inp, outp, n_threads, min_cont_len):
    if not os.path.exists('temp/Assembly_results'):
        os.mkdir('temp/Assembly_results')
    
    if os.path.exists('temp/Assembly_results/'+ outp):
        shutil.rmtree('temp/Assembly_results/'+ outp)
        
    outp = 'temp/Assembly_results/' + outp
    command = "./megahit -r {} --k-list 27,37,47,57,67,77,87 --min-contig-len {} -t {} -o {}".format(inp, min_cont_len, n_threads, outp)
    subprocess.check_call(command, shell=True)
    
def paired_assembly(inp1, inp2, outp, n_threads, min_cont_len):
    if not os.path.exists('temp/Assembly_results'):
        os.mkdir('temp/Assembly_results')
    
    if os.path.exists('temp/Assembly_results/'+ outp):
        shutil.rmtree('temp/Assembly_results/'+ outp)
        
    outp = 'temp/Assembly_results/' + outp
    command = "./megahit -1 {} -2 {} --k-list 27,37,47,57,67,77,87 --min-contig-len {} -t {} -o {}".format(inp1, inp2, min_cont_len, n_threads, outp)
    subprocess.check_call(command, shell=True)

def from_sra_to_cont(inp_name, out_name=None, QC=True, n_threads=1, min_cont_len=300):
    if out_name == None:
        out_name = inp_name
        if '/' in out_name:
            out_name = out_name.split('/')[-1]
        elif '//' in out_name:
            out_name = out_name.split('//')[-1]

    sra2fq(inp_name, out_name)
    fqs = os.listdir('temp/FQs/'+out_name)
    
    if QC == True:
        if not os.path.exists('temp/QC_reports/'+out_name):
            os.mkdir('temp/QC_reports/'+out_name)
        for fq in fqs:
            QC_report('temp/FQs/{}/{}'.format(out_name, fq), out_name + '/' + fq)
    
    if len(fqs) == 1:
        single_assembly('temp/FQs/{}/{}'.format(out_name, fqs[0]), str(out_name), n_threads=n_threads, min_cont_len=min_cont_len)
        
    elif len(fqs) == 2:
        paired_assembly('temp/FQs/{}/{}'.format(out_name, fqs[0]), 'temp/FQs/{}/{}'.format(out_name, fqs[1]), 
                        str(out_name), n_threads=n_threads, min_cont_len=min_cont_len)
            
            
def from_fq_to_cont(list_of_fq, out_name=None, QC=True, n_threads=1, min_cont_len=300):
    if out_name == None:
        out_name = list_of_fq[0]
        if '/' in out_name:
            out_name = out_name.split('/')[-1]
        elif '//' in out_name:
            out_name = out_name.split('//')[-1]

    if not os.path.exists('temp/FQs'):
        os.mkdir('temp/FQs')
    
    if not os.path.exists('temp/FQs/'+out_name):
        os.mkdir('temp/FQs/'+ out_name)
    
    for file in list_of_fq:
        file1 = file
        if '/' in file1:
            file1 = file1.split('/')[-1]
        elif '//' in file1:
            file1 = file1.split('//')[-1]
        shutil.copyfile(file, 'temp/FQs/{}/{}'.format(out_name, file1))

    fqs = os.listdir('temp/FQs/'+out_name)
    
    if QC == True:
        if not os.path.exists('temp/QC_reports/'+out_name):
            os.mkdir('temp/QC_reports/'+out_name)
        for fq in fqs:
            QC_report('temp/FQs/{}/{}'.format(out_name, fq), out_name + '/' + fq)
        
    if len(fqs) == 1:
        single_assembly('temp/FQs/{}/{}'.format(out_name, fqs[0]), str(out_name),n_threads=n_threads, min_cont_len=min_cont_len)
        
    elif len(fqs) == 2:
        paired_assembly('temp/FQs/{}/{}'.format(out_name, fqs[0]), 'temp/FQs/{}/{}'.format(out_name, fqs[1]), 
                        str(out_name), n_threads=n_threads, min_cont_len=min_cont_len)
    
    
def collect_all_res(target_name):
    if not os.path.exists('MeTarEnz_Results'):
        os.mkdir('MeTarEnz_Results')
    if not os.path.exists('MeTarEnz_Results/' + target_name):
        os.mkdir('MeTarEnz_Results/' + target_name)
    if os.path.exists('MeTarEnz_Results/' + target_name):
        shutil.rmtree('MeTarEnz_Results/' + target_name)
        os.mkdir('MeTarEnz_Results/' + target_name)

    shutil.move('temp/FQs/'+target_name, 'MeTarEnz_Results/' + target_name + '/' + 'FastQ_files')
    shutil.move('temp/Assembly_results/'+target_name, 'MeTarEnz_Results/' + target_name + '/' + 'assembly_results')
    
    if os.path.exists('temp/QC_reports/'+target_name):
        shutil.move('temp/QC_reports/'+target_name, 'MeTarEnz_Results/' + target_name + '/' + 'Quality_Control_results')
    
if __name__ == '__main__':
    from_sra_to_cont('SRR12198400', out_name='SRR')
    collect_all_res('SRR')
    from_sra_to_cont('DRR220391', out_name='DRR')
    collect_all_res('DRR')
    from_fq_to_cont(['p_file1_1.fq', 'p_file1_2.fq'], out_name='p_f')
    collect_all_res('p_f')
    # sra2fq('test_inputs/SRR12198400', 'SRR12198400')
    # fqs = os.listdir('FQs/SRR12198400')
    
    # for fq in fqs:
    #     QC_report('FQs/SRR12198400/'+fq, fq)
        
    # if len(fqs) == 1:
    #     outp = fqs[0][:fqs[0].find('_')]
    #     single_assembly('FQs/SRR12198400/' + fqs[0], outp, n_threads=1, min_cont_len=300)
        
    # elif len(fqs) == 2:
    #     outp = fqs[0][:fqs[0].find('_')]
    #     paired_assembly('FQs/SRR12198400/'+fqs[0], 'FQs/SRR12198400/'+fqs[1], outp, n_threads=1, min_cont_len=300)