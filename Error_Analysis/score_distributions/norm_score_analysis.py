import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import seaborn as sns
from scipy.stats import rankdata


def normalize(scores, thr, method="none"):
    """Restituisce (scores trasformati, threshold trasformato, etichetta asse x).
    La stessa trasformazione monotona viene applicata a scores e threshold."""
    s = np.asarray(scores, dtype=float)
    if method == "none":
        return s, thr, "scores"
    if method == "minmax":
        lo, hi = s.min(), s.max()
        return (s - lo) / (hi - lo), (thr - lo) / (hi - lo), "scores (min-max)"
    if method == "zthr":
        sd = s.std()
        return (s - thr) / sd, 0.0, "distanza dal cut-off (unità di std)"
    if method == "rank":
        r = rankdata(s) / len(s)
        thr_r = (s < thr).mean()  # percentile a cui cade il cut-off
        return r, thr_r, "percentile di rango"
    raise ValueError(f"norm sconosciuta: {method}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Plot the distribution of TP, TN, FP, FN')
    parser.add_argument('--data', required=True, type=str, help = 'Path to the input file containing predictions and scores')
    parser.add_argument('--out', type=str, required=True, help='Name of the output plot file')
    parser.add_argument('--threshold',type=float,required=True,help='Threshold used for the classification')
    parser.add_argument('--norm', choices=['none', 'minmax', 'zthr', 'rank'], default='none',
                        help='Normalizzazione dell\'asse x (applicata anche al threshold): '
                             'none | minmax [0,1] | zthr (distanza dal cut-off in std) | rank (percentile)')
    parser.add_argument('--method', type=str, default='VHM', help='Nome del metodo da mettere nel titolo')

    args=parser.parse_args()

    X = pd.read_csv(args.data, sep=',')

    pd.options.display.float_format = "{:.15f}".format



    X["category"] = np.select(
    [(X["true_class"] == 0) & (X["pred_class"] == 0),
     (X["true_class"] == 0) & (X["pred_class"] == 1),
     (X["true_class"] == 1) & (X["pred_class"] == 0),
     (X["true_class"] == 1) & (X["pred_class"] == 1)],
    ["TN", "FP", "FN", "TP"], default="NA")

    # normalizzazione dell'asse x (scores e threshold insieme)
    X["x"], threshold, xlabel = normalize(X["scores"], args.threshold, args.norm)
    output_file = args.out

    fig, ax = plt.subplots(figsize=(9, 5))


    sns.kdeplot(data=X, x="x", hue="category", hue_order=["TN", "FP", "FN", "TP"],
            common_norm=True, fill=True, alpha=0.3, log_scale=False, cut=0, ax=ax)
 

    ax.axvline(threshold, color="black", linestyle="--", linewidth=2)
    title = f"Distribution of {args.method} decision scores by confusion matrix category"
    if args.norm != "none":
        title += f"\n(x normalized: {args.norm})"
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_yscale("function", functions=(np.sqrt, np.square))
    #ax.set_xlim(-10, 20)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.close() 

    # for execute
    # python norm_score_analysis.py --data /home/um40/LB2_project_Group_1/Models/vhm/predictions.csv --out graph_vhm.png --threshold 9.268752841808581 --method VHM --norm zthr
    # python /home/um40/LB2_project_Group_1/Error_Analysis/norm_score_analysis.py --data/home/um40/LB2_project_Group_1/Models/SVM/SVM_pred_vs_test.csv --out graph_svm_norm.png --threshold 0 --method SVM --norm zthr