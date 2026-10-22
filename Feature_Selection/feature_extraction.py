import numpy as np
import matplotlib.pyplot as plt

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData 
from aaindex_parse import parse_aaindex_file
from collections import Counter
import argparse



#Amino acid scale: Transmembrane tendency ProtScale

dict_f['tt'] = {
    'A':  0.380,  # Ala
    'R': -2.570,  # Arg
    'N': -1.620,  # Asn
    'D': -3.270,  # Asp
    'C': -0.300,  # Cys
    'Q': -1.840,  # Gln
    'E': -2.900,  # Glu
    'G': -0.190,  # Gly
    'H': -1.440,  # His
    'I':  1.970,  # Ile
    'L':  1.820,  # Leu
    'K': -3.460,  # Lys
    'M':  1.400,  # Met
    'F':  1.980,  # Phe
    'P': -1.440,  # Pro
    'S': -0.530,  # Ser
    'T': -0.320,  # Thr
    'W':  1.530,  # Trp
    'Y':  0.490,  # Tyr
    'V':  1.460,  # Val
}

def profile_stat(profile):
    return max(profile),np.argmax(profile)/40



AA_keys="ACDEFGHIKLMNPQRSTVWY"
def count_aa(seq):

    dict_aa= Counter({aa: 0 for aa in AA_keys})

    for aa in seq:
        dict_aa[aa]+=1

    return list(dict_aa.values())



def profile(sequence, scale, w):
    # one value per residue; ends and unknown residues handled with the shrinking window
    L = len(sequence)
    d = w // 2
    profile = []
    for i in range(L):
        window = []
        for j in range(i - d, i + d + 1):
            if 0 <= j < L and sequence[j] in scale:     # inside the sequence AND a known residue
                window.append(scale[sequence[j]])
        profile.append(sum(window) / len(window) if window else np.nan) #append the mean 
    return np.array(profile)




def extract_f(seq,dict_f):

    #______kd_________
    kd=profile(seq[:40],dict_f['kd'],5)
    kd_max,kd_max_p=profile_stat(kd)
    kd_th=2 #decision from literature
    kd_len = 0

    if kd_max>=kd_th:

        i=int(kd_max_p*40)

        kd_len=0
        while i<len(kd) and kd[i]>=kd_th:
            kd_len+=1
            i+=1

        i=int(kd_max_p*40-1)
        while i>=0 and kd[i]>=kd_th:
            kd_len+=1
            i-=1

    kd_dekay=kd_max-kd[-1]

    #______charge_________
    charge_pos=profile(seq[:40],dict_f['p'],3)
    max_charge_pos,max_charge_pos_p=profile_stat(charge_pos)

    charge_neg=profile(seq[:40],dict_f['n'],3)
    max_charge_neg,max_charge_neg_p=profile_stat(charge_neg)

    #_____nfh_&_nfb____________
    nfh=profile(seq[:40],dict_f['nfh'],9)
    max_nfh,max_nfh_p=profile_stat(nfh)

    nfb=profile(seq[:40],dict_f['nfb'],9)
    max_nfb,max_nfb_p=profile_stat(nfb)

    #____res_vol_____
    res_vol_mean=np.mean(profile(seq[:40],dict_f['v'],5))

    #____signal_seq_h_potential
    signal=profile(seq[:40],dict_f['s'],9)
    signal_max,signal_max_p=profile_stat(signal)

    #___flexibility___
    flexibility_mean=np.mean(profile(seq[:40],dict_f['f'],9))


    #___mut____
    mut_mean=np.mean(profile(seq[:40],dict_f['rm'],9))


    #___N-ter_comp_____
    n_ter_comp=count_aa(seq)

    return [kd_max,kd_max_p,kd_len,kd_dekay,max_charge_pos,max_charge_pos_p,max_charge_neg,max_charge_neg_p,max_nfh,max_nfh_p,max_nfb,max_nfb_p,res_vol_mean,signal_max,signal_max_p,flexibility_mean,mut_mean]+n_ter_comp



if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Feature Selection")
    parser.add_argument('--scale_selected_file',type=str,default='/home/um40/LB2_project_Group_1/Feature_Selection/AAindex_scales.txt')
    parser.add_argument('--fasta',type=str,default='/home/um40/LB2_project_Group_1/von_Heijne_method/merged.fasta')
    parser.add_argument('--aaindex_file',type=str,default='/home/um40/LB2_project_Group_1/Feature_Selection/aaindex1',help='download the file from aaindex.com')
    args=parser.parse_args()

    with open(args.scale_selected_file,'r') as f:
        dict_f={}
        for line in f:
            line=line.strip().split()
            dict_f[line[1]]=parse_aaindex_file(args.aaindex_file,line[0])


    with open('feature.csv','w') as feature_file:
        print('id,max_kd,p_max_kd,kd_len,kd_dakay,max_charge_pos,max_charge_pos_p,max_charge_neg,max_charge_neg_p,nfh_max,nfh_p,nfb_max,nfb_p,res_vol_mean,signal_max,signal_max_p,flexibility_mean,mut_mean,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y',file=feature_file)
    
        with open(args.fasta,'r') as fasta:
            for line in fasta:
                line=line.strip()
                if line.startswith('>'):
                    id=line[1:]
                    seq=fasta.readline().strip()
                    f=[id]+extract_f(seq,dict_f)
                    print(*f,sep=',',file=feature_file)
