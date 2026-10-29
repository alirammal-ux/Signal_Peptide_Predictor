import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

#AAAAAAAAAAAAAAAAAAAAAAAAAAA

#    /\_/\          MIAOOOOOO
#   ( o.o )   <-- Il gatto della qualità del codice veglia su di te.
#    > ^ <

#FIXME: better to be pass from command line
#HERE = Path(__file__).resolve().parent   # Feature_Selection/
#ROOT = HERE.parent                       # repository root


argparser = argparse.ArgumentParser(description='This script is used to analyze the predictions of the model on the benchmark set')


def performance_metrics(y_test,y_pred):
    '''This function is used to store the performance metrics of the model on the benchmark set'''
    mcc= matthews_corrcoef(y_test, y_pred)
    precision= precision_score(y_test, y_pred)
    recall= recall_score(y_test, y_pred)
    f1=f1_score(y_test, y_pred)
    cm=confusion_matrix(y_test, y_pred)
    bal_acc=balanced_accuracy_score(y_test, y_pred)

    return mcc, f1, cm, precision, recall, bal_acc


def labeling_function(df):
    TN=df.loc[(df['true_class']==0) & (df['pred_class']==0)]['UniprotID'].tolist()

    TP=df.loc[(df['true_class']==1) & (df['pred_class']==1)]['UniprotID'].tolist()

    FN=df.loc[(df['true_class']==1) & (df['pred_class']==0)]['UniprotID'].tolist()

    FP=df.loc[(df['true_class']==0) & (df['pred_class']==1)]['UniprotID'].tolist()

    return TN, TP, FN, FP


def plot_taxa_distribution(FP,FN,data):
    FP_taxa=data.loc[data['ID'].isin(FP)]['kingdom'].value_counts()
    FN_taxa=data.loc[data['ID'].isin(FN)]['kingdom'].value_counts()

    taxa_df = pd.DataFrame({'False Positives': FP_taxa, 'False Negatives': FN_taxa})
    taxa_df_normalized = taxa_df.div(taxa_df.sum(axis=1), axis=0) * 100 

    taxa_df_normalized.plot(kind='bar', color=['red', 'blue'], alpha=0.5)


    plt.title('Distribution of False Positives and False Negatives across Taxa', weight='bold')
    plt.xlabel('Taxa')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig('taxa_distribution.png', dpi=300)
    plt.show()



def Helix_presence(FP,data):
    FP_helix = data.loc[data['ID'].isin(FP) & (data['early_trans'] == True)]['early_trans'].value_counts()
    FP_No_helix = data.loc[data['ID'].isin(FP) & (data['early_trans'] == False)]['early_trans'].value_counts()

    counts = pd.DataFrame({'Helix Presence': FP_helix, 'No Helix Presence': FP_No_helix})

    counts.plot(kind='bar', color=['red', 'orange'], alpha=0.5)

    plt.title('Distribution of False Positives across Helix Presence', weight='bold')
    plt.xlabel('Helix Presence')
    plt.ylabel('Count')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig('helix_distribution.png', dpi=300)
    plt.show()



if __name__=="__main__":
    argparser.add_argument('--output_model',  type=str, required=True, help='Path to the output model file')
    argparser.add_argument('--original_df', type=str, required=True, help='Path to the original dataframe file containg the information from uniprot')
    args = argparser.parse_args()
    df=pd.read_csv(args.output_model)
    data=pd.read_csv(args.original_df)



