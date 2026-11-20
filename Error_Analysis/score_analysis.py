

import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import seaborn as sns


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Plot the distribution of TP, TN, FP, FN')
    parser.add_argument('--data', required=True, type=str, help = 'Path to the input file containing predictions and scores')
    parser.add_argument('--out', type=str, help='Name of the output plot file')
    parser.add_argument('--threshold',type=float,help='Threshold used for the classification')
    
    args=parser.parse_args()

    X = pd.read_csv(args.data, sep=',')

    pd.options.display.float_format = "{:.15f}".format
    print(X.head())



    X["category"] = np.select(
    [(X["true_class"] == 0) & (X["pred_class"] == 0),
     (X["true_class"] == 0) & (X["pred_class"] == 1),
     (X["true_class"] == 1) & (X["pred_class"] == 0),
     (X["true_class"] == 1) & (X["pred_class"] == 1)],
    ["TN", "FP", "FN", "TP"])

    threshold = args.threshold
    output_file = args.out

    fig, ax = plt.subplots(figsize=(9, 5))


    sns.kdeplot(data=X, x="scores", hue="category", hue_order=["TN", "FP", "FN", "TP"],
            common_norm=False, fill=True, alpha=0.3, log_scale=False, cut=0, ax=ax)
 

    ax.axvline(threshold, color="black", linestyle="--", linewidth=2)
    ax.set_title("Distribution of VHM method decision scores by confusion matrix category")
    ax.set_yscale("function", functions=(np.sqrt, np.square))
    #ax.set_xlim(-10, 20)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close() 

    # for execute
    # python score_analysis.py --data /home/um40/LB2_project_Group_1/Models/vhm/predictions.csv --out graph_vhm.png --threshold 9.268752841808581
    # python score_analysis.py --data /home/um40/LB2_project_Group_1/Models/SVM/SVM_pred_vs_test.csv --out graph_svm.png --threshold 0






