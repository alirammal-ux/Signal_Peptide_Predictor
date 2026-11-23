import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve,confusion_matrix,matthews_corrcoef, f1_score,precision_score, recall_score,accuracy_score,ConfusionMatrixDisplay
import argparse
from pathlib import Path
import matplotlib.pyplot as plt 


AA_keys="ACDEFGHIKLMNPQRSTVWY"
N_TER=90
HERE = Path(__file__).resolve().parent   # von_Heijne_method/
ROOT = HERE.parent.parent                # root del repository


def count(fasta,data,j=None,final=False):
    set1=list()
    set2=list()


    with open(fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                UniprotID = line.rstrip()[1:]

                #take the fold 
                fold = data.loc[data['ID'] == UniprotID, 'fold'].values[0] # estraggo la riga corrispondente ad UniprotID e ricavo il valore di fold per determinare a che set appartiene

                #take the class
                clas=data.loc[data['ID'] == UniprotID, 'class'].values[0]

                #take the context
                cleavage_pos = data.loc[data['ID'] == UniprotID, 'cleavage_pos'].values[0]


                start_pos = cleavage_pos - 13
                end_pos = cleavage_pos + 2

                seq = fasta.readline().strip()


                if final:
                    if fold==-1:
                        set2.append((seq,clas,UniprotID))
                    else:
                        if clas==1: 
                            set1.append(seq[start_pos:end_pos])

                else:

                    #lista per split 
                    fold_train=[1,2,3,4,5]
                    fold_train.remove(j)

                    if fold in fold_train and clas==1: 
                        set1.append(seq[start_pos:end_pos])
                    elif fold==j:
                        set2.append((seq,clas))


    return set1,set2

#pesudocount

def create_pswm(train_set,freq):
    pspm=np.ones([15,20])

    for seq in train_set:
        for i in range(len(seq)):
            pos_aa=AA_keys.find(seq[i])
            # residui non standard (X, U, ...): find() restituisce -1, che li conterebbe come Y
            if pos_aa == -1:
                continue
            pspm[i,pos_aa]+=1

    # senza residui saltati ogni riga somma a N+20: equivale a dividere per (len(train_set)+20)
    pspm/=pspm.sum(axis=1, keepdims=True)


    return np.log2(pspm/freq)


def score_seq(seq,pswm):

    best_win_score,pos_best=float('-inf'),0
    final_pos=min(len(seq),N_TER)-15
    
    #scorro sulla val_seq
    for i in range(0,final_pos+1):
        score=0
        
        #scorro sulla pswm
        for row in range(0,len(pswm[:,0])):


            #trovo col aa nella pswm 
            col=AA_keys.find(seq[row+i])

            # residuo non standard (X, U, ...): non contribuisce allo score
            if col == -1:
                continue
            
            #aggiungo allo score della window, il valore dell aa nella pswm
            score+=pswm[row,col]

        #esaurisco la window e aggiorno il miglior score solo se è il migliore 
        if best_win_score < score:
            best_win_score=score
            pos_best=i

    return best_win_score,pos_best


def opt_th(y_val , y_val_scores):
    # y: 0=NO_SP, 1=SP ; scores from PSWM scanning on validation data
    precision, recall, thresholds = precision_recall_curve(y_val , y_val_scores)
    f1 = (2 * precision * recall) / (precision + recall)

    # se precision+recall=0 f1 è NaN: nanargmax lo ignora, argmax restituirebbe proprio il NaN
    optimal_threshold = thresholds[np.nanargmax(f1)]

    return optimal_threshold



def metrics_test(th,y_test_scores,y_test):
    y_pred = [int(s >= th) for s in y_test_scores]

    cm= confusion_matrix(y_test, y_pred)
    mcc=matthews_corrcoef(y_test,y_pred)
    f1=f1_score(y_test,y_pred)
    precision=precision_score(y_test,y_pred)
    recall=recall_score(y_test,y_pred)
    acc=accuracy_score(y_test,y_pred)

    return (mcc,f1,precision,recall,acc),cm




if __name__=='__main__':
    
    #argparse
    parser = argparse.ArgumentParser(description="Calcolo metriche su test set")
    parser.add_argument('--df_back',type=str,default=str(ROOT/'Data_analysis'/'aa_composition_plots'/'background_freq.tsv'))
    parser.add_argument('--fasta',type=str,default=str(ROOT/'von_Heijne_method/merged.fasta'))
    parser.add_argument('--data',type=str,default=str(ROOT/'Data_preparation'/'Cleaned_merged_dataset.tsv'),help='Cleaned_merged_dataset')
    parser.add_argument('--out_dir',type=str,default=str(HERE))
    parser.add_argument('--cm',type=str,default='n',help='y if you want a cm for run,default=n')
    args=parser.parse_args()
    out_dir=Path(args.out_dir)
    df_background=pd.read_table(args.df_back)
    df_background.sort_values(by='aa',inplace=True)
    fasta=args.fasta
    df=pd.read_table(args.data)

    freq=df_background['freq_back'].to_numpy()

    list_th=[]
    


    df_metrics = pd.DataFrame(columns=['mcc', 'f1', 'precision', 'recall', 'acc','tn', 'fp', 'fn', 'tp', 'threshold'])
    
    #######VALiDATION#####################
    for run in range(1,6):
        
        train_set,val_set=count(fasta,df,run)
        pswm=create_pswm(train_set,freq)

        y_val_scores=[]
        y_val=[]
        
        for (seq,clas) in val_set:
            y_val_scores.append(score_seq(seq,pswm)[0])
            y_val.append(clas)


        th=opt_th(y_val , y_val_scores)

        list_th.append(th)

    #######TEST#####################
    train_set,test_set_with_ID=count(fasta,df,final=True)
    pswm=create_pswm(train_set,freq)
    y_test_scores=[]
    y_test=[]

    #mean of fold_threshold
    best_threshold=np.mean(list_th)
            
    for (seq,clas,_) in test_set_with_ID:
        y_test_scores.append(score_seq(seq,pswm)[0])
        y_test.append(clas)
    

    m,cm=metrics_test(best_threshold,y_test_scores,y_test)
    tn, fp, fn, tp = cm.ravel()
    mcc,f1,precision,recall,acc=m
    df_metrics.loc[len(df_metrics)]=mcc,f1,precision,recall,acc,tn, fp, fn, tp, best_threshold

    df_pred = pd.DataFrame({
        'UniprotID': [item[2] for item in test_set_with_ID],
        'true_class': [item[1] for item in test_set_with_ID],
        'pred_class': [int(score >= best_threshold) for score in y_test_scores],
        'scores':y_test_scores
    })
    df_pred.to_csv(out_dir / f'predictions.csv', index=False)

    if args.cm == 'y':
        disp = ConfusionMatrixDisplay(cm, display_labels=["non-SP", "SP"])
        fig, ax = plt.subplots(figsize=(5, 4))
        disp.plot(ax=ax, cmap="Blues", colorbar=False, values_format="d")
        ax.set_title(f"Confusion matrix – final")
        fig.savefig(out_dir / f"confusion_matrix_final.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


    #save df 
    df_metrics.to_csv(out_dir/'vom_heijne_metrics.csv',index=False)
