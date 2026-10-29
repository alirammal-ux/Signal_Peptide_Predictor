import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import argparse
from collections import Counter
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

argparser = argparse.ArgumentParser(description="Analyze the results of the SVM test set.")

##### AA COMPOSITION at N-TERM ######
N_TERM = 40
AA_keys="ACDEFGHIKLMNPQRSTVWY"

groups = [
    "nonpolar",   # A
    "polar",      # C
    "negative",   # D
    "negative",   # E
    "aromatic",   # F
    "nonpolar",   # G
    "positive",   # H
    "nonpolar",   # I
    "positive",   # K
    "nonpolar",   # L
    "nonpolar",   # M
    "polar",      # N
    "nonpolar",   # P
    "polar",      # Q
    "positive",   # R
    "polar",      # S
    "polar",      # T
    "nonpolar",   # V
    "aromatic",   # W
    "aromatic",   # Y
]

def count_aa(seq, dict_aa): #FIXME: better import from the first file in which is defined (be careful because this version is more general than the original one => return dict without normalization)
    for aa in seq:
        if aa in AA_keys: #skip unknown residues
            dict_aa[aa]+=1
    return dict_aa


def compute_aa_composition(seqs):
    aa_dict= Counter({aa: 0 for aa in AA_keys}) # initialization of the dictionary

    for seq in seqs:
        aa_dict=count_aa(seq,aa_dict) #update the dictionary for each FN seq beloging to benchmark set

    # Normalize
    for aa in aa_dict.keys():
        aa_dict[aa]/=(N_TERM*len(seqs)) #FIXME: in this way we are normalizing considering also the unknown residues
    
    return aa_dict.values()


def generate_comp_bargraph(FN_comp, TP_comp, TP_comp_train, output_dir):
    # Collect all required data in a df
    df_comp = {
        'aa': list(AA_keys),
        'FN_test': FN_comp,
        'TP_test': TP_comp,
        'TP_train': TP_comp_train,
        'group': groups
    }
    df_comp = pd.DataFrame(df_comp)

    # 2. Order the df with respect group
    df_comp.sort_values(by='group', inplace=True)

    #sns.set_theme(style="white", font_scale=1.05)
    #plt.rcParams['font.family'] = 'DejaVu Sans'

    palette = {
        'FN_test': "#2E5EAA", 
        'TP_test': "#D64545", 
        'TP_train': "#4CAF7A"
    }

    # 3. Melt the dataframe
    # transform the 3 columns specified in values_vars (= FN_comp, TP_comp, TP_comp_train) 
    # into a single column (= Frequency) where the corresponding labels are mapped in a new column (= metric, which takes values from {FN_comp, TP_comp, TP_comp_train})
    # => from 3 columns of numeric values we obtain just two columns, one of numbers and the other of labels (=FN_comp, TP_comp, TP_comp_train)
    df_melted = df_comp.melt(
        id_vars=['aa', 'group'], 
        value_vars=['FN_test', 'TP_test', 'TP_train'],
        var_name='metric', 
        value_name='Frequency'
    )

    # Ordine degli amminoacidi raggruppati
    aa_order = df_comp.sort_values('group')['aa'].tolist()
    aa_to_group = dict(zip(df_comp['aa'], df_comp['group']))

    fig, ax = plt.subplots(figsize=(20, 7.5))

    # 4. Passaggio di hue='metric' anziché 'aa'
    sns.barplot(
        data=df_melted, 
        x='aa', 
        y='Frequency', 
        hue='metric',
        palette=palette, 
        edgecolor='white', 
        linewidth=1.2, 
        ax=ax,
        order=aa_order
    )

    # Etichette dentro/sopra le barre (corretti i segni meno mancanti)
    for patch in ax.patches:
        height = patch.get_height()
        if pd.isna(height) or height == 0:
            continue
        ax.text(
            patch.get_x() + patch.get_width() / 2, # Posizione X al centro della barra
            #height - 0.5 if height > 0.5 else height + 0.002,
            height / 2, # Posizione Y a metà altezza della barra 
            f'{height:.3f}', # Formattazione a 2 cifre decimali
            ha='center', va='center',
            rotation=90, # Ruotato di 90 gradi (puoi metterlo a 0 se lo preferisci orizzontale)
            fontsize=12, color='white', fontweight='bold'
        )

    # Pulizia estetica
    ax.set_title('AA composition of FN_test vs TP_test vs TP_train', fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel('Residue Type', fontsize=12, labelpad=10)
    ax.set_ylabel('Frequency', fontsize=12, labelpad=10)

    sns.despine(left=True, bottom=True)
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    ax.set_axisbelow(True)

    ax.legend(title='', frameon=False, loc='upper right', fontsize=11)

    # --- Graffe di gruppo sotto l'asse x (corretti i segni meno mancanti) ---
    def add_group_brackets(ax, aa_order, aa_to_group, y_line=-0.14, y_text=-0.20,
                            tick_height=0.025, color='black', linewidth=2.8):
        trans = ax.get_xaxis_transform()
        groups_seq = [aa_to_group[a] for a in aa_order]
        segments = []
        start = 0
        for i in range(1, len(groups_seq) + 1):
            if i == len(groups_seq) or groups_seq[i] != groups_seq[start]:
                segments.append((start, i - 1, groups_seq[start]))  # <--- Corretto i - 1
                start = i

        for start_idx, end_idx, group in segments:
            x_start = start_idx - 0.4  # <--- Corretto start_idx - 0.4
            x_end = end_idx + 0.4

            ax.plot([x_start, x_end], [y_line, y_line],
                    transform=trans, color=color, linewidth=linewidth,
                    clip_on=False, solid_capstyle='butt')
            ax.plot([x_start, x_start], [y_line, y_line + tick_height],
                    transform=trans, color=color, linewidth=linewidth, clip_on=False)
            ax.plot([x_end, x_end], [y_line, y_line + tick_height],
                    transform=trans, color=color, linewidth=linewidth, clip_on=False)

            x_center = (x_start + x_end) / 2
            ax.text(x_center, y_text, group, transform=trans,
            ha='center', va='top', fontsize=14, color=color,
            clip_on=False)

    add_group_brackets(ax, aa_order, aa_to_group, color='black')

    plt.tight_layout()
    #plt.show()
    plt.savefig(f'{output_dir}/Bargraph_aa_composition.png')


#### OTHER FEATURES DISTRIBUTION ####
feature_names_dict = {
    "max_kd": "Maximum Kyte-Doolittle Hydrophobicity",
    "p_max_kd": "Position of Maximum Kyte-Doolittle Hydrophobicity",
    "kd_len": "Length of Hydrophobic Region",
    "kd_dakay": "Kyte-Doolittle Hydrophobicity decay",
    "max_charge_pos": "Maximum Positive Charge",
    "max_charge_pos_p": "Position of Maximum Positive Charge",
    "max_charge_neg": "Maximum Negative Charge",
    "max_charge_neg_p": "Position of Maximum Negative Charge",
    "nfh_max": "Maximum NFH Score",
    "nfh_p": "Position of Maximum NFH Score",
    "nfb_max": "Maximum NFB Score",
    "nfb_p": "Position of Maximum NFB Score",
    "res_vol_mean": "Mean Residue Volume",
    "signal_max": "Maximum Signal Peptide Score",
    "signal_max_p": "Position of Maximum Signal Peptide Score",
    "flexibility_mean": "Mean Protein Flexibility",
    "mut_mean": "Mean Mutation Score",
    "A": "Alanine Frequency",
    "C": "Cysteine Frequency",
    "D": "Aspartic Acid Frequency",
    "E": "Glutamic Acid Frequency",
    "F": "Phenylalanine Frequency",
    "G": "Glycine Frequency",
    "H": "Histidine Frequency",
    "I": "Isoleucine Frequency",
    "K": "Lysine Frequency",
    "L": "Leucine Frequency",
    "M": "Methionine Frequency",
    "N": "Asparagine Frequency",
    "P": "Proline Frequency",
    "Q": "Glutamine Frequency",
    "R": "Arginine Frequency",
    "S": "Serine Frequency",
    "T": "Threonine Frequency",
    "V": "Valine Frequency",
    "W": "Tryptophan Frequency",
    "Y": "Tyrosine Frequency"
}

def feature_distributions_boxplot(FN_id_test, TP_id_test, TP_id_train, df_features, feature, output_dir):
    FN_test = df_features.loc[df_features["id"].isin(FN_id_test)][feature].tolist() 
    TP_test = df_features.loc[df_features["id"].isin(TP_id_test)][feature].tolist()
    TP_train = df_features.loc[df_features["id"].isin(TP_id_train)][feature].tolist()

    plt.figure(figsize=(8, 5))
    plt.boxplot([FN_test, TP_test, TP_train], labels=['FN test', 'TP test', 'TP train'])

    plt.title(f'Comparison of {feature_names_dict[feature]} distributions')
    plt.ylabel(feature_names_dict[feature])
    plt.grid(True, linestyle='--', alpha=0.6)
    #plt.show()
    plt.savefig(f'{output_dir}/Boxplot_{feature}_dist.png')







if __name__ == '__main__':
    argparser.add_argument('--output_model',  type=str, required=True, help='Path to the output model file')
    argparser.add_argument('--fasta', type=str, required=True, help='Path to the fasta file containing the sequences of the von Heijne test set')
    argparser.add_argument('--output_dir', type=str, required=True, help='path to the output directory')
    argparser.add_argument('--original_df', type=str, required=True, help='Path to the original dataframe')
    argparser.add_argument('--features_csv', type=str, required=True, help='Path to the csv file containing all the features')
    args = argparser.parse_args()

    df_labels = pd.read_csv(args.output_model)
    data = pd.read_csv(args.original_df)
    df_features = pd.read_csv(args.features_csv)


    TN_idx = df_labels.loc[(df_labels['true_class']==0) & (df_labels['pred_class']==0)]['UniprotID'].tolist()
    TP_idx = df_labels.loc[(df_labels['true_class']==1) & (df_labels['pred_class']==1)]['UniprotID'].tolist()
    FN_idx = df_labels.loc[(df_labels['true_class']==1) & (df_labels['pred_class']==0)]['UniprotID'].tolist()
    FP_idx = df_labels.loc[(df_labels['true_class']==0) & (df_labels['pred_class']==1)]['UniprotID'].tolist()

    TP_idx_train = data.loc[(data["fold"]!=-1)&(data["class"]==1)]['ID'].tolist()

    #### AA COMPOSITION at N-TERM ####

    FN_seq=[] # False negative sequences (till pos 40=k) beloging to bechmark/test set 
    TP_seq=[] # True positive sequences (till pos 40=k) beloging to bechmark/test set 
    TP_seq_train=[] #True positive sequences (till pos 40=k) beloging to training set

    with open(args.fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                UniprotID = line.strip()[1:]
                seq=fasta_file.readline().strip()
                if UniprotID in FN_idx:
                    FN_seq.append(seq[0:N_TERM])
                elif UniprotID in TP_idx:
                    TP_seq.append(seq[0:N_TERM])
                elif UniprotID in TP_idx_train:
                    TP_seq_train.append(seq[0:N_TERM])

    FN_comp = compute_aa_composition(FN_seq)
    TP_comp = compute_aa_composition(TP_seq)
    TP_comp_train = compute_aa_composition(TP_seq_train)

    generate_comp_bargraph(FN_comp, TP_comp, TP_comp_train, args.output_dir)

    #### SIGNAL PEPTIDE LENGTH DISTRIBUTION ####

    sp_length_FN = data.loc[data["ID"].isin(FN_idx)]['cleavage_pos'].tolist()
    sp_length_TP = data.loc[data["ID"].isin(TP_idx)]['cleavage_pos'].tolist()
    sp_length_TP_train = data.loc[data["ID"].isin(TP_idx_train)]['cleavage_pos'].tolist()

    # Linegraph of the three distribution by using KDE (Kernel Density Estimation)
    # kde computes the esimate kernel density, and this is necessary for correctly comparing the three lists in the same graph
    # because they have different length (= diff #seq)
    plt.figure(figsize=(10, 6))

    sns.kdeplot(sp_length_FN, label='FN (Test Set)', palette="rocket")
    sns.kdeplot(sp_length_TP, label='TP (Test Set)', palette="rocket")
    sns.kdeplot(sp_length_TP_train, label='TP (Train Set)', palette="rocket")

    plt.title('Comparison of signal peptide length distributions', fontsize=14)
    plt.xlabel('Length (aa)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f'{args.output_dir}/Linegraph_sp_length_dist.png')

    # Boxplots of the three distributions
    plt.figure(figsize=(8, 5))
    plt.boxplot([sp_length_FN, sp_length_TP, sp_length_TP_train], labels=['sp_length_FN', 'sp_length_TP', 'sp_length_TP_train'])

    plt.title('Comparison of signal peptide length distributions')
    plt.ylabel('Peptide length')
    plt.grid(True, linestyle='--', alpha=0.6)
    #plt.show()
    plt.savefig(f'{args.output_dir}/Boxplots_sp_length_dist.png')


    #### OTHER FEATURES DISTRIBUTIONS ####

    feature_distributions_boxplot(FN_idx, TP_idx, TP_idx_train, df_features, 'max_kd', args.output_dir)
    feature_distributions_boxplot(FN_idx, TP_idx, TP_idx_train, df_features, 'signal_max', args.output_dir)    











