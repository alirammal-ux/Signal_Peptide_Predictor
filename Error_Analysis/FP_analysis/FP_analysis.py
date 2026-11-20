import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path


HERE = Path(__file__).resolve().parent  # FP_analysis/
ROOT = HERE.parent.parent               # repository root
N_SPECIES = 4

def performance_FPR_analysis(y_true, y_pred, TP, TN, FP, FN, data, model, output_dir):

    ##### PERFORMACE EVALUATION #####
    mcc= matthews_corrcoef(y_true, y_pred)
    precision= precision_score(y_true, y_pred)
    recall= recall_score(y_true, y_pred)
    f1=f1_score(y_true, y_pred)
    bal_acc=balanced_accuracy_score(y_true, y_pred)
    cm=confusion_matrix(y_true, y_pred)

    # Generate heatmap for confusion matrix
    # 4. Generazione della heatmap
    plt.figure(figsize=(5, 4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['Non-SP', 'SP'], yticklabels=['Non-SP', 'SP'], cbar=False) # --> add normalize='all' to obtain proportions instead of count values
    plt.title(f'Confusion matrix ({model})')
    plt.ylabel('True labels')
    plt.xlabel('Predicted labels')
    plt.savefig(f'{output_dir}/{model}_cm.png')
    #plt.show()

    ##### FPR ANALYSIS #####
    # 1. Extract the uniprotID corresponding to TN, TP, FN, FP wich are transmembrane (TM) proteins from both 
    TN_tm = data.loc[data['ID'].isin(TN) & data['early_trans']==True]['ID'].tolist()
    FP_tm = data.loc[data['ID'].isin(FP) & data['early_trans']==True]['ID'].tolist()
    TP_tm = data.loc[data['ID'].isin(TP) & data['early_trans']==True]['ID'].tolist()
    FN_tm = data.loc[data['ID'].isin(FN) & data['early_trans']==True]['ID'].tolist()
    # print(f'|TN_tm|={len(TN_tm)}\n|FP_tm|={len(FP_tm)}\n|TP_tm|={len(TP_tm)}\n|FN_tm|={len(FN_tm)}') 
    # The fact that no transmembran helix is present in the positive dataset (TP+FN) is ok
    
    # 1. Compute FPR (both tot and for only TM)
    FPR_tm=len(FP_tm)/(len(FP_tm)+len(TN_tm))
    FPR_tot=len(FP)/(len(TP)+len(TN))
    

    # Write result in output file
    with open(f'{output_dir}/{model}_performance_&_FPR_analysis.txt', 'w') as txt_file:
        print(f'PERFORMANCE EVALUATION\nmcc: {mcc}\nprecision: {precision}\nrecall: {recall}\nf1: {f1}\nbalanced accuracy: {bal_acc}\nconfusion matrix: {cm}', file=txt_file)
        print(f'\nFALSE POSITIVE RATE ANALYSIS\nFPR_tot: {FPR_tot}\nFPR_tm: {FPR_tm}', file=txt_file)

    return mcc, f1, precision, recall, bal_acc



def pieplot_norm_distribution(FP, FN, TP, TN, feature, type, model, output_dir):
    # 1. extract the number of entries for each feature value (ex: feature=kingdom/organism_name) in the benchmark/test fold
    FP_grouped = data.loc[data['ID'].isin(FP)][feature].value_counts()
    FN_grouped = data.loc[data['ID'].isin(FN)][feature].value_counts()
    TP_grouped = data.loc[data['ID'].isin(TP)][feature].value_counts()
    TN_grouped = data.loc[data['ID'].isin(TN)][feature].value_counts()  

    # 2. generate pieplot
    fig, ax= plt.subplots(2,2, figsize=(20,16), facecolor = 'white')

    ax[0, 0].pie(FP_grouped.values, labels=FP_grouped.index, autopct='%1.1f%%',  textprops={'fontsize': 16})
    ax[0, 0].set_title(f"False Positive (|FP|={len(FP)})", fontsize=20)

    ax[0, 1].pie(TN_grouped.values, labels=TN_grouped.index, autopct='%1.1f%%', textprops={'fontsize': 16})
    ax[0, 1].set_title(f"True Negative (|TN|={len(TN)})", fontsize=20)

    ax[1, 0].pie(FN_grouped.values, labels=FN_grouped.index, autopct='%1.1f%%', textprops={'fontsize': 16})
    ax[1, 0].set_title(f"False Negative (|FN|={len(FN)})", fontsize=20)

    ax[1, 1].pie(TP_grouped.values, labels=TP_grouped.index, autopct='%1.1f%%',  textprops={'fontsize': 16})
    ax[1, 1].set_title(f"True Positive (|TP|={len(TP)})", fontsize=20)

    fig.suptitle(f'Distribution of False Positives and False Negatives across {type} ({model})', weight='bold', fontsize=20, y=0.95)
    
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(f'{output_dir}/{model}_pieplot_{type}_dist.png')
    #plt.show()


def spider_plot(df, metrics, output_path, models, title="Benchmarking of models' performances"):
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
                marker="o", markersize=4, label=models[i])


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

    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    argparser = argparse.ArgumentParser(description='This script is used to analyze the predictions of the model on the benchmark set')

    argparser.add_argument('--output_model',  nargs='+', default=[ROOT/'Models'/'SVM'/'SVM_pred_vs_test.csv', ROOT/'Models'/'vhm'/'predictions.csv', ROOT/'Models'/'CNN'/'test_predictions.csv'], help='Path to the output model file')
    argparser.add_argument('--original_df', type=str, default=ROOT/'Data_preparation'/'Cleaned_merged_dataset.tsv', help='Path to the original dataframe file containg the information from uniprot')
    argparser.add_argument('--output_dir', type=str, default=HERE, help='Path to the output directory where the plots will be saved')
    argparser.add_argument('--model_name', nargs='+', default=['svm', 'vhm', 'cnn'], help='Name of the model used')

    args = argparser.parse_args()

    output_model_files = args.output_model # --> default: [svm_file, vhm_file]
    data = pd.read_csv(args.original_df, sep='\t')
    output_dir = args.output_dir
    models_names = args.model_name
    df_metrics = pd.DataFrame(columns=['mcc', 'f1', 'precision', 'recall', 'balanced acc'])

    # sanity check
    assert len(output_model_files) == len(models_names), "Error: the number of output_model files does not correspond to the number of model_name inserted"

    N = len(models_names) # --> number of models to be analyzed 
    for i in range(0,N):
        df = pd.read_csv(output_model_files[i])
        model = models_names[i]

        ####################################
        #### PERFORMANCE & FPR ANALYSIS ####
        ####################################

        # 1.1) Extract the UniprotID corresponding to TN, TP; FN; FP in the benchmark/test setof the IDs 
        TN = df.loc[(df['true_class']==0) & (df['pred_class']==0)]['UniprotID'].tolist()
        TP = df.loc[(df['true_class']==1) & (df['pred_class']==1)]['UniprotID'].tolist()
        FN = df.loc[(df['true_class']==1) & (df['pred_class']==0)]['UniprotID'].tolist()
        FP = df.loc[(df['true_class']==0) & (df['pred_class']==1)]['UniprotID'].tolist()

        # 1.2) Compute performance metrics, FPR tot, FPR only for TM protein
        mcc, f1, precision, recall, bal_acc  = performance_FPR_analysis(y_true=df['true_class'], y_pred=df['pred_class'], TP=TP, TN=TN, FP=FP, FN=FN, data=data, model=model, output_dir=output_dir) 
        df_metrics.loc[len(df_metrics)] = mcc, f1, precision, recall, bal_acc

    
        ###################################
        ##### FP DISTRIBUTION in TAXA #####
        ###################################

        # 2) Generate the four pieplot (FP vs TN, FN vs TP)
        pieplot_norm_distribution(FP=FP, FN=FN, TP=TP, TN=TN, feature='kingdom', type='Taxa', model=model, output_dir=output_dir)  

        ###################################
        ##### FP DISTRIBUTION in TAXA #####
        ###################################

        # 3.1) Replace longest species name with a shorter version
        data['organism_name'] = data['organism_name'].replace('Saccharomyces cerevisiae (strain ATCC 204508 / S288c)', 'Saccharomyces cerevisiae')
        data['organism_name'] = data['organism_name'].replace('Schizosaccharomyces pombe (strain 972 / ATCC 24843)', 'S. pombe')

        # 3.2) Find the N species most frequents and rename all the other as "other"
        top = data['organism_name'].value_counts().head(N_SPECIES).index
        data.loc[~data['organism_name'].isin(top), 'organism_name'] = 'Other'
        # In "df[~df", the tilde operator means NOT!, and so is used in this scenario to filter out the top n organisms from the dataframe. 
        
        # 3.3) Generate the four pieplot (FP vs TN, FN vs TP)
        pieplot_norm_distribution(FP=FP, FN=FN, TP=TP, TN=TN, feature='organism_name', type='Species', model=model, output_dir=output_dir) 
    
    # Generate the spider plot to visualize the benchmarking
    spider_plot(df_metrics, ["mcc", "f1", "precision", "recall", "balanced acc"], models=models_names, output_path=f'{output_dir}/benchmarking_spider_plot.png')
            

    



