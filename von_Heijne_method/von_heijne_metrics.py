import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve,confusion_matrix,matthews_corrcoef, f1_score,precision_score, recall_score,accuracy_score
import argparse
import os


AA_keys="ACDEFGHIKLMNPQRSTVWY" 


def count(fasta,data,j):
    cleavage_site_train=list()
    seq_val=list()
    seq_test=list()


    with open(fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                UniprotID = line.rstrip()[1:]

                #take the fold 
                fold = data.loc[data['ID'] == UniprotID, 'fold'].values[0] -1 # estraggo la riga corrispondente ad UniprotID e ricavo il valore di fold per determinare a che set appartiene

                #take the class
                clas=data.loc[data['ID'] == UniprotID, 'class'].values[0]

                #take the context
                cleavage_pos = data.loc[data['ID'] == UniprotID, 'cleavage_pos'].values[0]


                start_pos = cleavage_pos - 13
                end_pos = cleavage_pos + 2

                seq = fasta.readline().strip()

                #lista per split 
                fold_train=[(i+j)%5 for i in range(2,5)]
                

                if fold in fold_train and clas==1: 
                    cleavage_site_train.append(seq[start_pos:end_pos])
                elif fold==(j+1)%5:
                    seq_val.append((seq,clas))
                elif fold==j:
                    # print('è entrato')
                    seq_test.append((seq,clas))

                    
    return cleavage_site_train,seq_val,seq_test


#pesudocount

def create_pswm(train_set,freq):
    pspm=np.ones([15,20])

    for seq in train_set:
        for i in range(len(seq)):
            pos_aa=AA_keys.find(seq[i])
            pspm[i,pos_aa]+=1

    pspm/=(len(train_set)+20)


    return np.log2(pspm/freq)


def score_seq(seq,pswm):

    best_win_score,pos_best=float('-inf'),0
    final_pos=min(len(seq),90)-15
    
    #scorro sulla val_seq
    for i in range(0,final_pos):
        score=0
        
        #scorro sulla pswm
        for row in range(0,len(pswm[:,0])):


            #trovo col aa nella pswm 
            col=AA_keys.find(seq[row+i])
            
            #aggiungo allo score della window, il valore dell aa nella pswm
            score+=pswm[row,col]

        #esaurisco la window e aggiorno il miglior score solo se è il migliore 
        if best_win_score < score:
            best_win_score=score
            pos_best=i

    return best_win_score,pos_best


def opt_th(y_val , y_val_scores,y_test_scores):
    # y: 0=NO_SP, 1=SP ; scores from PSWM scanning on validation data
    precision, recall, thresholds = precision_recall_curve(y_val , y_val_scores)
    f1 = (2 * precision * recall) / (precision + recall)

    optimal_threshold = thresholds[np.argmax(f1)]
    # classify the test set with the chosen threshold
    # y_pred = [int(s >= optimal_threshold) for s in y_test_scores]

    return optimal_threshold

def metrics_test(th,y_test_scores,y_test):
    y_pred = [int(s >= th) for s in y_test_scores]
    cm=confusion_matrix(y_test,y_pred)

    mcc=matthews_corrcoef(y_test,y_pred)
    f1=f1_score(y_test,y_pred)
    precision=precision_score(y_test,y_pred)
    recall=recall_score(y_test,y_pred)
    acc=accuracy_score(y_test,y_pred)

    return mcc,f1,precision,recall,acc,cm





if __name__=='__main__':

    
    #argparse
    parser = argparse.ArgumentParser(description="Calcolo metriche su test set")
    parser.add_argument('--df_back',type=str,default='/home/um40/LB2_project_Group_1/Data_analysis/aa_composition_plots/background_freq.tsv')
    parser.add_argument('--fasta',type=str,default='/home/um40/LB2_project_Group_1/von_Heijne_method/merged.fasta')
    parser.add_argument('--data',type=str,default='/home/um40/LB2_project_Group_1/Data_preparation/Cleaned_merged_dataset.tsv')
    parser.add_argument('--out',type=str,default='/home/um40/LB2_project_Group_1/von_Heijne_method')
    args=parser.parse_args()

    dir_path=os.path(args.out)
    df_background=pd.read_table(args.df_back)
    fasta=args.fasta
    df=pd.read_table(args.data)

    freq=df_background['freq_back'].to_numpy()

    df_metrics = pd.DataFrame(columns=['mcc', 'f1', 'precision', 'recall', 'acc', 'cm'])

    for run in range(0,5):
        train_set,val_set,test_set=count(fasta,df,run)
        pswm=create_pswm(train_set,freq)

        y_val_scores=[]
        y_val=[]
        for (seq,clas) in val_set:
            y_val_scores.append(score_seq(seq,pswm)[0])
            y_val.append(clas)

        y_test_scores=[]
        y_test=[]
        for (seq,clas) in test_set:
            y_test_scores.append(score_seq(seq,pswm)[0])
            y_test.append(clas)


        th=opt_th(y_val , y_val_scores,y_test_scores)
        df_metrics.loc[len(df_metrics)]=metrics_test(th,y_test_scores,y_test)


    #salvo df 
    df_metrics.to_csv(f'{dir_path}vom_heijne_metrics.csv')

    se_matrics=df_metrics[df_metrics.columns[:-1]].std(axis=0)/np.sqrt(len(df_metrics))
    avg_metrics=df_metrics[df_metrics.columns[:-1]].mean(axis=0)

    
    for metrics in df_metrics.columns[:-1]:
        print(f'{metrics} = {avg_metrics[metrics]} ± {se_matrics[metrics]}')
