import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve,confusion_matrix,matthews_corrcoef, f1_score,precision_score, recall_score,accuracy_score,ConfusionMatrixDisplay
import argparse
from pathlib import Path
import matplotlib.pyplot as plt 


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

    cm= confusion_matrix(y_test, y_pred,normalize='true')
    mcc=matthews_corrcoef(y_test,y_pred)
    f1=f1_score(y_test,y_pred)
    precision=precision_score(y_test,y_pred)
    recall=recall_score(y_test,y_pred)
    acc=accuracy_score(y_test,y_pred)

    return (mcc,f1,precision,recall,acc),cm


def spider_plot(df, metrics, out_path, title="Performance across folds"):
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 7), subplot_kw=dict(polar=True))

    # Sfondo bianco
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Primo asse in alto, senso orario
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    colors = plt.cm.Set1.colors   # palette satura: rosso, blu, verde, viola, arancio...

    for i, row in df[metrics].iterrows():
        vals = row.tolist() + [row.iloc[0]]
        ax.plot(angles, vals, color=colors[i % len(colors)], linewidth=1.8, alpha=0.9,
                marker="o", markersize=4, label=f"Fold {i}")

    mean = df[metrics].mean().tolist()
    mean += mean[:1]
    ax.plot(angles, mean, color="black", linewidth=2.5, label="Mean")
    ax.fill(angles, mean, color="black", alpha=0.1)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics)
    ax.set_ylim(0, 1)
    ax.spines["polar"].set_color("#cccccc")
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=8, color="#666666")
    ax.yaxis.grid(True, color="#999999", linestyle="--", linewidth=0.8)   # cerchi concentrici
    ax.xaxis.grid(True, color="#999999", linestyle="-", linewidth=0.8)    # raggi
    ax.set_axisbelow(True)                                                # griglia sotto le linee dei dati
    ax.spines["polar"].set_color("#999999")
    ax.tick_params(axis="x", pad=10)
    ax.set_title(title, pad=25)

    # Legenda sotto, centrata
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=9, frameon=False)

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)



if __name__=='__main__':

    
    #argparse
    parser = argparse.ArgumentParser(description="Calcolo metriche su test set")
    parser.add_argument('--df_back',type=str,default='/home/um40/LB2_project_Group_1/Data_analysis/aa_composition_plots/background_freq.tsv')
    parser.add_argument('--fasta',type=str,default='/home/um40/LB2_project_Group_1/von_Heijne_method/merged.fasta')
    parser.add_argument('--data',type=str,default='/home/um40/LB2_project_Group_1/Data_preparation/Cleaned_merged_dataset.tsv',help='Cleaned_merged_dataset')
    parser.add_argument('--out_dir',type=str,default='/home/um40/LB2_project_Group_1/von_Heijne_method')
    parser.add_argument('--cm',type=str,default='n',help='y if you want a cm for run,default=n')
    parser.add_argument('--spider',type=str,default='n',help='y if you want a spider for run,default=n')
    args=parser.parse_args()

    out_dir=Path(args.out_dir)
    df_background=pd.read_table(args.df_back)
    fasta=args.fasta
    df=pd.read_table(args.data)

    freq=df_background['freq_back'].to_numpy()

    df_metrics = pd.DataFrame(columns=['mcc', 'f1', 'precision', 'recall', 'acc','tn', 'fp', 'fn', 'tp'])

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
        m,cm=metrics_test(th,y_test_scores,y_test)
        tn, fp, fn, tp = cm.ravel()
        mcc,f1,precision,recall,acc=m
        df_metrics.loc[len(df_metrics)]=mcc,f1,precision,recall,acc,tn, fp, fn, tp

        if args.cm == 'y':
            disp = ConfusionMatrixDisplay(cm, display_labels=["non-SP", "SP"])
            fig, ax = plt.subplots(figsize=(5, 4))
            disp.plot(ax=ax, cmap="Blues", colorbar=False, values_format=".2%")
            ax.set_title(f"Confusion matrix – run {run}")
            fig.savefig(out_dir / f"confusion_matrix{run}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            



    #save df 
    df_metrics.to_csv(out_dir/'vom_heijne_metrics.csv')

    se_metrics=df_metrics[df_metrics.columns].std(axis=0)/np.sqrt(len(df_metrics))
    avg_metrics=df_metrics[df_metrics.columns].mean(axis=0)

    if args.spider=='y':
        spider_plot(df_metrics, ["mcc", "f1", "precision", "recall", "acc"], out_dir / "spider_plot.png")

    with open(out_dir / 'k_fold_stats.csv', 'w') as file:
        print('metric,mean,se', file=file)
        for metric in df_metrics.columns:
            print(f'{metric},{avg_metrics[metric]:.4f},{se_metrics[metric]:.4f}', file=file)

