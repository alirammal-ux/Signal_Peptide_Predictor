import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import argparse
from pathlib import Path


#HERE = Path(__file__).resolve().parent   # Feature_Selection/
#ROOT = HERE.parent                       # repository root


argparser = argparse.ArgumentParser(description='This script is used to analyze the predictions of the model on the benchmark set')




#AAAAAAAAAAAAAAAAAAAAAAAAAAA

#    /\_/\          MIAOOOOOO
#   ( o.o )   <-- Il gatto della qualità del codice veglia su di te.
#    > ^ <
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
plt.rcParams["figure.figsize"] = (15,9)
## extraction of the dataframes:
#output of the model df

if __name__=="__main__":
    argparser.add_argument('--output_model',  type=str, required=True, help='Path to the output model file')
    argparser.add_argument('--original_df', type=str, required=True, help='Path to the original dataframe file containg the information from uniprot')
    argparser.add_argument('--output_dir', type=str, required=True, help='Path to the output directory where the plots will be saved')
    argparser.add_argument('--model_name', type=str, required=True, help='Name of the model used')

    args = argparser.parse_args()
    df=pd.read_csv(args.output_model)
    data=pd.read_csv(args.original_df, sep='\t')
    output_dir=args.output_dir
    model_name=args.model_name



# performance extraction

# EVALUATION OF THE PERFORMANCE:

def performance_metrics(y_test,y_pred):
    '''This function is used to store the performance metrics of the model on the benchmark set'''
    mcc= matthews_corrcoef(y_test, y_pred)
    precision= precision_score(y_test, y_pred)
    recall= recall_score(y_test, y_pred)
    f1=f1_score(y_test, y_pred)
    cm=confusion_matrix(y_test, y_pred)
    bal_acc=balanced_accuracy_score(y_test, y_pred)

    return mcc, f1, cm, precision, recall, bal_acc
performance=performance_metrics(df['true_class'], df['pred_class'])
performance
#### extraction of the ids
#In the following cell we extracted the IDs of the elements of the confusion matrix, in order to be able to do the analysis. 
def labeling_function(df):
    TN=df.loc[(df['true_class']==0) & (df['pred_class']==0)]['UniprotID'].tolist()

    TP=df.loc[(df['true_class']==1) & (df['pred_class']==1)]['UniprotID'].tolist()

    FN=df.loc[(df['true_class']==1) & (df['pred_class']==0)]['UniprotID'].tolist()

    FP=df.loc[(df['true_class']==0) & (df['pred_class']==1)]['UniprotID'].tolist()

    return TN, TP, FN, FP

TN, TP, FN, FP = labeling_function(df)
tn_trans=data.loc[data['ID'].isin(TN) & data['early_trans']==True]['ID'].tolist()
print(len(tn_trans))
fp_trans=data.loc[data['ID'].isin(FP) & data['early_trans']==True]['ID'].tolist()

print(len(fp_trans))

tp_trans=data.loc[data['ID'].isin(TP) & data['early_trans']==True]['ID'].tolist()
print(len(tp_trans))
fn_trans=data.loc[data['ID'].isin(FN) & data['early_trans']==True]['ID'].tolist()
print(len(fn_trans))

## sono meno perche sono solo sul test set, mentre i primi sono sul dataset completo. 

#Already from here the fact that no transmembran elix is present in the positive dataset  (TP+FN) is ok.
data['early_trans'].value_counts()
TN
fpr_tm=len(fp_trans)/(len(fp_trans)+len(tn_trans))
fpr_tot=len(FP)/(len(TP)+len(TN))
print(f'fpr_tot: {fpr_tot}\nfpr_tm: {fpr_tm}')
#The FPR of the transmembrane proteins is 10 times higher than the FPR of the total dataset. This was expected since the residues of protein with a transmembrane secondary structure at N terminal, is more prone to be confused by the model as the hydrophobic region of the signal peptide. 
## Distribution analysis:
data.kingdom.unique()
#extraction of number of entries for taxa in the test set
fungi_test=len(data[(data['kingdom']=='Fungi') & (data['fold']==-1)]['ID'].to_list())
metazoa_test=len(data[(data['kingdom']=='Metazoa') & (data['fold']==-1)]['ID'].to_list())
plants_test=len(data[(data['kingdom']=='Plants') & (data['fold']==-1)]['ID'].to_list())
other_test=len(data[(data['kingdom']=='other') & (data['fold']==-1)]['ID'].to_list())

taxa=data.loc[data["fold"] ==-1].groupby('kingdom').size()
print(fungi_test, metazoa_test, plants_test, other_test)
print(taxa)
FP_taxa=data.loc[data['ID'].isin(FP)]['kingdom'].value_counts()
FN_taxa=data.loc[data['ID'].isin(FN)]['kingdom'].value_counts()

def plot_taxa_distribution(FP,FN,data,output_dir, model_name):
    FP_taxa=data.loc[data['ID'].isin(FP)]['kingdom'].value_counts()
    FN_taxa=data.loc[data['ID'].isin(FN)]['kingdom'].value_counts()

    taxa_df = pd.DataFrame({'False Positives': FP_taxa, 'False Negatives': FN_taxa})
    taxa_df_normalized = taxa_df.div(taxa_df.sum(axis=1), axis=0) * 100 

    ax= taxa_df.plot(kind='bar', color=['red', 'blue'], alpha=0.5)

    for p in ax.patches:
        h = p.get_height()
        if int(h)==0:
            pass
        else:
            ax.annotate(f'{int(h)}%',
                    (p.get_x() + p.get_width() / 2, h / 2),
                    ha='center', va='center')

    plt.title('Distribution of False Positives and False Negatives across Taxa', weight='bold')
    plt.xlabel('Taxa')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{model_name}_taxa_distribution.png', dpi=900, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.show()

def plot_taxa_norm_distribution(FP,FN,data,output_dir, model_name):
    FP_taxa=data.loc[data['ID'].isin(FP)]['kingdom'].value_counts()
    FN_taxa=data.loc[data['ID'].isin(FN)]['kingdom'].value_counts()

    print(FN_taxa)
    print(taxa)

    FP_taxa = FP_taxa / taxa.sum() * 100
    FN_taxa = FN_taxa / taxa.sum() * 100

    taxa_df = pd.DataFrame({'False Positives': FP_taxa, 'False Negatives': FN_taxa})
    #taxa_df_normalized = taxa_df.div(taxa_df.sum(axis=1), axis=0) * 100 
    taxa_df_normalized = taxa_df

    ax=taxa_df_normalized.plot(kind='bar', color=['red', 'blue'], alpha=0.5)

    for p in ax.patches:
            h = p.get_height()
            if int(h)==0:
                pass
            else:
                ax.annotate(f'{int(h)}%',
                        (p.get_x() + p.get_width() / 2, h / 2),
                        ha='center', va='center')


    plt.title('Distribution of False Positives and False Negatives across Taxa', weight='bold')
    plt.xlabel('Taxa')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{model_name}_taxa_normalized_distribution.png', dpi=900, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.show()



def Helix_presence(FP,data,output_dir, model_name):
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
    plt.savefig(f'{output_dir}/{model_name}_helix_distribution.png', dpi=900, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.show()



def Helix_presence_norm(FP,data,output_dir, model_name):
    FP_helix = data.loc[data['ID'].isin(FP) & (data['early_trans'] == True)]['early_trans'].value_counts()
    FP_No_helix = data.loc[data['ID'].isin(FP) & (data['early_trans'] == False)]['early_trans'].value_counts()

    Number_helix = data.loc[data['early_trans'] == True]['early_trans'].value_counts()
    Number_No_helix = data.loc[data['early_trans'] == False]['early_trans'].value_counts()

    FP_helix = FP_helix / Number_helix * 100
    FP_No_helix = FP_No_helix / Number_No_helix * 100

    counts = pd.DataFrame({'Helix Presence': FP_helix, 'No Helix Presence': FP_No_helix})

    ax= counts.plot(kind='bar', color=['red', 'orange'], alpha=0.5)
    for p in ax.patches:
        h = p.get_height()
        if h != 0:
            ax.annotate(
                f'{float(h):.2f}%',
                (p.get_x() + p.get_width() / 2, h),
                ha='center',
                va='bottom',
                xytext=(0, 3),
                textcoords='offset points'
            )

    plt.title('Distribution of False Positives across Helix Presence', weight='bold')
    plt.xlabel('Helix Presence')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{model_name}_helix_distribution_norm.png', dpi=900)
    plt.show()



### TAXA:
#We can observe that the distribution of the false positive and negatives in all the taxas are pretty much similar, changing at most of 2%. 


#plot_taxa_norm_distribution(FP,FN,data)
def plot_taxa_all_distribution(FP,FN, TP, TN, data,output_dir, model_name):
    Correctly_predicted = data.loc[data['ID'].isin(TP or TN)]['kingdom'].value_counts()
    FP_taxa=data.loc[data['ID'].isin(FP)]['kingdom'].value_counts()
    FN_taxa=data.loc[data['ID'].isin(FN)]['kingdom'].value_counts()

    

    taxa_df = pd.DataFrame({'False Positives': FP_taxa, 'False Negatives': FN_taxa, 'Correctly Predicted': Correctly_predicted})

    taxa_df_normalized = taxa_df.div(taxa_df.sum(axis=1), axis=0) * 100

    ax= taxa_df_normalized.plot(kind='bar', color=['red', 'blue', 'green'], alpha=0.5)

    for p in ax.patches:
        h = p.get_height()
        if int(h)==0:
            pass
        else:
            ax.annotate(f'{int(h)}%',
                    (p.get_x() + p.get_width() / 2, h / 2),
                    ha='center', va='center')


 
    plt.title('Distribution of False Positives and False Negatives across Kingdoms', weight='bold')
    plt.xlabel('Kingdoms')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.figure(figsize=(20, 12))
    plt.savefig(f'{output_dir}/{model_name}_kingdom_all_distribution.png', dpi=900)
    plt.show()

def plot_species_all_distribution(FP,FN, TP, TN, data,output_dir, model_name):

    plt.rcParams["figure.figsize"] = (15,9)


    n = 3
    top = data['organism_name'].value_counts().head(n).index

    # In "df[~df", the tilde operator means NOT!, and so is used in this scenario to filter out the top n organisms from the dataframe.  
    data.loc[~data['organism_name'].isin(top), 'organism_name'] = 'Other'
    data['organism_name'] = data['organism_name'].replace('Saccharomyces cerevisiae (strain ATCC 204508 / S288c)', 
        'Saccharomyces cerevisiae'
    )
    data['organism_name'] = data['organism_name'].replace('Schizosaccharomyces pombe (strain 972 / ATCC 24843)', 
        'S. pombe'
    )

    Correctly_predicted = data.loc[data['ID'].isin(TP or TN)]['organism_name'].value_counts()
    FP_taxa=data.loc[data['ID'].isin(FP)]['organism_name'].value_counts()
    FN_taxa=data.loc[data['ID'].isin(FN)]['organism_name'].value_counts()

    print(FN_taxa)
    print(FP_taxa)
    print(Correctly_predicted)


    

    taxa_df = pd.DataFrame({'False Positives': FP_taxa, 'False Negatives': FN_taxa, 'Correctly Predicted': Correctly_predicted})

    taxa_df_normalized = taxa_df.div(taxa_df.sum(axis=1), axis=0) * 100

    ax= taxa_df_normalized.plot(kind='bar', color=['red', 'blue', 'green'], alpha=0.5)

    for p in ax.patches:
        h = p.get_height()
        if int(h)==0:
            pass
        else:
            ax.annotate(f'{int(h)}%',
                    (p.get_x() + p.get_width() / 2, h / 2),
                    ha='center', va='center')


 
    plt.title('Distribution of False Positives and False Negatives across Species', weight='bold')
    plt.xlabel('Species')
    plt.ylabel('percentage')
    plt.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{model_name}_species_all_distribution.png', dpi=900, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.show()

plot_species_all_distribution(FP,FN, TP, TN, data,output_dir, model_name)
plot_taxa_all_distribution(FP,FN, TP, TN, data,output_dir, model_name)
### HELIX:
Helix_presence(FP,data,output_dir, model_name)
Helix_presence_norm(FP,data,output_dir, model_name)



