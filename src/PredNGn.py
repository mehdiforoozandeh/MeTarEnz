from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet
import numpy as np
from featgen import featgen
import joblib


def feature_extraction(seq):
    f_vec = featgen.Extract_all(seq)
    return f_vec


def read_fasta(file_name):
    from Bio import SeqIO
    seqsdic = {}
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    for fasta in fasta_sequences:
        seqsdic[fasta.id] = str(fasta.seq)
    return seqsdic


def obtain_seq_from_ent(ent):
    from Bio import Entrez, SeqIO
    Entrez.email = "mehdiforoozandehsh@gmail.com"

    try:
        handle = Entrez.efetch(db="protein", id=ent, rettype="gp")
        record = SeqIO.read(handle, "gb")
        seq = str(record.seq)

    except:
        from pydpi import pypro
        seq = str(pypro.GetProteinSequence(ent))

    if not _verify_alphabet(Seq(seq, IUPAC.protein)):
        print('Inserted Entry Is Not Valid!')
    return seq


def seq_repair(seq):
    forbidden_chars = ['B','J','O','U','X','Z','@','_','!','#','$','%','^','&','*',
                       '(',')','<','>','?','/','|','}','{','~',':','-', '\n', ' ']
    for f in forbidden_chars:
            while f in seq:
                seq = seq.replace(f,'')
    return seq


def ph_f_selection(f_vec):
    with open('Models/pHdupidx.txt', 'r') as conf:
        idx = conf.readlines()[0].split(',')
    for i in range(len(idx)):
        idx[i] = int(idx[i])
    idx = np.array(idx)
    f_vec = f_vec[:, np.sort(idx)]
    scaler = joblib.load('Models/pH_Scaler.sav')
    f_vec = scaler.transform(f_vec)
    selector = joblib.load('Models/pH_Selector.sav')
    x = selector.transform(f_vec)
    pca = joblib.load('Models/pH_PCA.sav')
    x = pca.transform(x)
    return x


def temp_f_selection(f_vec):
    with open('Models/tempdupidx.txt', 'r') as conf:
        idx = conf.readlines()[0].split(',')
    for i in range(len(idx)):
        idx[i] = int(idx[i])
    idx = np.array(idx)
    f_vec = f_vec[:, np.sort(idx)]
    scaler = joblib.load('Models/temp_Scaler.sav')
    f_vec = scaler.transform(f_vec)
    selector = joblib.load('Models/temp_Selector.sav')
    x = selector.transform(f_vec)
    pca = joblib.load('Models/temp_PCA.sav')
    x = pca.transform(x)
    return x


def phpredict(x):
    model = joblib.load('Models/pH_model.sav')
    outp = model.predict(x).reshape(1, -1)

    return outp


def tempredict(x):
    model = joblib.load('Models/temp_model.sav')
    outp = model.predict(x).reshape(1, -1)

    return outp


# _____TOOL1: Prediction_____
def single_prediction(seq):
    if len(seq) > 20 and _verify_alphabet(Seq(seq, IUPAC.protein)):
        seq = seq_repair(seq)
        full_features = np.array(feature_extraction(seq))
        ph_x = ph_f_selection(full_features)
        opt_ph = phpredict(ph_x)[0]
        temp_x = temp_f_selection(full_features)
        opt_temp = tempredict(temp_x)[0]
        return opt_temp[0], opt_ph[0]
    
    else:
        return 'Not a Valid Protein Sequence!', 'Not a Valid Protein Sequence!'


def fasta_prediction(file_name):
    seqsdic = read_fasta(file_name)
    result = []
    for seq_id, seq in seqsdic.items():
        seq = seq_repair(seq)
        full_features = np.array(feature_extraction(seq))
        ph_x = ph_f_selection(full_features)
        opt_ph = phpredict(ph_x)[0][0]
        temp_x = temp_f_selection(full_features)
        opt_temp = tempredict(temp_x)[0][0]
        result.append([str(seq_id), opt_ph, opt_temp])
        
    return result