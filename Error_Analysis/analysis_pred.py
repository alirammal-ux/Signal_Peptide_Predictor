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

if __name__=="__main__":
    argparser.add_argument('--output_model',  type=str, required=True, help='Path to the output model file')
    args = argparser.parse_args()
    df=pd.read_csv(args.output_model)


