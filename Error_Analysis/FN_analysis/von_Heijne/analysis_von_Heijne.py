import pandas as pd
import argparse

argparser = argparse.ArgumentParser(description="Analyze the results of the von Heijne test set.")





def analyze_von_heijne_FN(fasta_file, prediction, data, output_dir):
    cleavage_site_FN=list()
    cleavage_site_positive=list()

    data = pd.read_csv(data, sep='\t')
    pred = pd.read_csv(prediction)

    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                UniprotID = line.rstrip()[1:]
                FN = pred.loc[(pred['UniprotID'] == UniprotID) & (pred['true_class'] == 1) & (pred['pred_class'] == 0)]['UniprotID']
                positive_training = data.loc[(data['ID'] == UniprotID) & (data['class'] == 1)]['ID']
                cleavage_pos = data.loc[data['ID'] == UniprotID, 'cleavage_pos'].values[0]
                start_pos = cleavage_pos - 13
                end_pos = cleavage_pos + 2
                seq = fasta.readline()





                if FN.any():
                    cleavage_site_FN.append(seq[start_pos:end_pos])

                elif positive_training.any():
                    cleavage_site_positive.append(seq[start_pos:end_pos])

                else:
                    continue

                
    with open(f'{output_dir}/cleavage_region_positive.txt', 'w') as txt:
        for seq in cleavage_site_positive:
            print(seq, file=txt)

    with open(f'{output_dir}/cleavage_region_FN.txt', 'w') as txt:
        for seq in cleavage_site_FN:
            print(seq, file=txt)

    
if __name__=="__main__":
    argparser.add_argument('--output_model',  type=str, required=True, help='Path to the output model file')
    argparser.add_argument('--fasta', type=str, required=True, help='Path to the fasta file containing the sequences of the von Heijne test set')
    argparser.add_argument('--output_dir', type=str, required=True, help='path to the output directory')
    argparser.add_argument('--original_df', type=str, required=True, help='Path to the original dataframe containing the cleavage positions')
    args = argparser.parse_args()
    pred=pd.read_csv(args.output_model)
    data=pd.read_csv(args.original_df)

    analyze_von_heijne_FN(args.fasta, args.output_model, args.original_df, args.output_dir)

