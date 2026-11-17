'''In this file we import the features and train the SVM model using the features selected by the Random Forest importance. This is done with a 5 fold cv, then the model will be tested on the benchmark  set'''

import time
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import matthews_corrcoef,f1_score, confusion_matrix, precision_score, recall_score, balanced_accuracy_score
from pathlib import Path
import argparse
import joblib

RANDOM_STATE = 0
HERE = Path(__file__).resolve().parent   # Feature_Selection/
ROOT = HERE.parent.parent                    # repository root


######### DEFINITION OF THE FUNCTION TO TRAIN AND VALIDATE THE MODEL #######


def train_and_test(X_train, y_train, X_test, y_test, new_model, c, gamma):
    '''If this function is called by choose_setting, the X_test and Y_test are the validation set. Otherwise this function can be called using the benchmark for the actual tesing of the model'''
    # train a classifier on (X_train, y_train) and return its mcc on (X_test, y_test)
    
    if X_train.shape[1] == 0:
        # no feature selected: the best we can do is predict the most frequent class
        # (its MCC is 0: returning the accuracy would not be comparable with the other settings)
        majority = np.bincount(y_train).argmax()
        y_pred = np.full(len(y_test), majority)
        return matthews_corrcoef(y_test, y_pred)
    
    # 1. learn the scaling (mean, std of each column) on the TRAINING data only
    scaler = StandardScaler().fit(X_train) # --> zeta scores = subtract the mean and divide by std (dependency from data, we compute mean & std on training deta => NOT USE TEST DATA TO COMPUTE THESE METRICS)
    # 2. train the classifier on the scaled training data
    model = new_model(c,gamma)
    model.fit(scaler.transform(X_train), y_train)
    # 3. scale the test data with the SAME mean and std, predict, measure mcc
    y_pred = model.predict(scaler.transform(X_test))
    

    return y_pred

def performance_metrics(y_test,y_pred):
    '''This function is used to store the performance metrics of the model on the benchmark set'''
    mcc= matthews_corrcoef(y_test, y_pred)
    precision= precision_score(y_test, y_pred)
    recall= recall_score(y_test, y_pred)
    f1=f1_score(y_test, y_pred)
    cm=confusion_matrix(y_test, y_pred)
    bal_acc=balanced_accuracy_score(y_test, y_pred)

    return mcc, f1, cm, precision, recall, bal_acc
    


######## CROSS VALIDATION FUNCTION#####
def choose_setting(select, grid, X, y, new_model,test_fold,fold):
    '''Takes as input the selection function, grid of hyperpameters, data and label, the model and the test fold. If it is not nested the test fold must be the benchmark fold (in our case -1).
    It returns the best setting and the mean validation mcc of every setting.'''

    # inner cross-validation on (X_train, y_train) only; returns the best setting 
    # and the mean validation mcc of every setting
    train_fold=[k for k in range(1,6) if k!=test_fold] 
    mean_val_mcc = {}
    for setting in grid: # --> for each set of hyperparam values in the grid
        k, c, gamma = setting
        val_mccs = []
        for j in train_fold:

            in_val=fold[fold==j].index.to_list()
            in_train = fold[(fold != j) & (fold != test_fold) & (fold > 0)].index.to_list()
    
            # selection sees the inner-training part only
            cols = select(X.iloc[in_train], y.iloc[in_train], k) #--> select the best features
             
            y_pred = train_and_test(X.iloc[in_train, cols], y.iloc[in_train],
                                 X.iloc[in_val, cols],   y.iloc[in_val], new_model,c,gamma) #--> train a new model with the subset of best features, then test it on validation fold
            mcc=performance_metrics(y.iloc[in_val], y_pred)[0]
            val_mccs.append(mcc)
        mean_val_mcc[setting] = np.mean(val_mccs)

    best = max(grid, key=mean_val_mcc.get)   # ties -> the first setting in grid (= smallest k)
    # if we have two combination of features with same mcc, always choose the one with fewer features (= less model's params, less overfit)
    return best, mean_val_mcc


def select_rf(X_train, y_train, k):
    rf = RandomForestClassifier(n_estimators=100, random_state=RANDOM_STATE)
    rf.fit(X_train, y_train)
    ranking = np.argsort(rf.feature_importances_)[::-1]   # most important first (= gini importance)
    return ranking[:k]


def new_svm(c,gamma):
    return SVC(kernel="rbf", C=c, gamma=gamma)


def model_eval(name, family, select, grid, new_model,X,y,fold,verbose=True, store=True):
    benchmark_fold=-1
    start = time.time()

    test_idx=fold[fold==benchmark_fold].index.to_list()
    train_idx=fold[(fold!=benchmark_fold) &( fold>0)].index.to_list()

    X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
    X_test,  y_test  = X.iloc[test_idx],  y.iloc[test_idx]

    #this function is used to choose the best setting on the training set, then train the model on the training set and test it on the benchmark set
    best, curve = choose_setting(select, grid, X, y, new_model,benchmark_fold,fold) #--> choose best setting (= #best_features to keep) on TRAIN SET

    k, c, gamma = best

    # STEP 2 - redo the selection on the whole outer-training part
        
    cols = select(X_train, y_train, k) #--> given the #best_features (just a number!), select which are the #best most significant features on TRAIN TEST

    scaler = StandardScaler().fit(X_train.iloc[:, cols])
    model = new_model(c, gamma)
    model.fit(scaler.transform(X_train.iloc[:, cols]), y_train)


    # STEP 3 - train on the outer-training part, test ONCE on the outer-test fold
    # y_pred = train_and_test(X_train.iloc[:, cols], y_train, X_test.iloc[:, cols], y_test, new_model, c, gamma) #--> train a new model with the subset of best features, then test it on validation fold
    y_pred = model.predict(scaler.transform(X_test.iloc[:, cols]))


    return y_test,y_pred,test_idx,k,c,gamma,model,scaler,cols

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='final training and evaluation of the svm')
    parser.add_argument('--X',default=ROOT/'SVM_method'/'feature.csv',type=str)
    parser.add_argument('--data',default=ROOT/'Data_preparation'/'Cleaned_merged_dataset.tsv',type=str)
    parser.add_argument('--out',default=HERE/'SVM_pred_vs_test.csv',type=str)
    parser.add_argument('--v',default='y',type=str,help='create a file SVM_performance.txt in the currend dir with a verbose description of the model performance')
    args=parser.parse_args()
    


    X=pd.read_csv(args.X) 

    # sort both tables by protein ID and reset the index, so that position i is the same protein in X, y and fold
    # (the fold indices are used as positions with .iloc in choose_setting and nested_cv)
    X = X.sort_values(by='id').reset_index(drop=True) #--> the values are sorted by the protein ID, so that the position i is the same protein in X, y and fold

    df=pd.read_csv(args.data,sep='\t')  

    df = df.sort_values(by='ID').reset_index(drop=True)

    ## we need to extract from df only the class and fold columns, and from X only the features (drop the id column)
    y= df['class']
    X= X.drop(columns='id', axis=1) #eliminate the id column from X, because it is not a feature, only used to synncronize with df
    fold= df['fold']

    feature_names=X.columns.to_list()



    ######### DEFINITION OF THE GRID OF PARAMETERS TO TEST #######

    #FEATURES K (number of features to select)
    K_GRID=[i for i in range(11,len(feature_names),2)]

    # C VALUES
    C_GRID=[1,4,10,16]

    #GAMMA VALUES dependent on k
    GRID=[]
    for k in K_GRID:
        for c in C_GRID:
            GRID.append((k,c,np.round(0.5/k,4))) 
            GRID.append((k,c,np.round(1/k,4))) 
            GRID.append((k,c,np.round(2/k,4))) 
            GRID.append((k,c,np.round(4/k,4))) 

    

    y_test, y_pred, test_idx, k, c, gamma, final_model, final_scaler, features = model_eval("SVM", "SVM", select_rf, GRID, new_svm,X,y,fold,verbose=True, store=True)

    mcc, f1, cm, precision, recall, bal_acc= performance_metrics(y_test, y_pred)

    # --- save model into pkl file  ---
    model_output_path = HERE / 'trained_svm_model.pkl'
    scaler_output_path = HERE / 'trained_scaler.pkl'
    
    joblib.dump(final_model, model_output_path)
    joblib.dump(final_scaler, scaler_output_path)
    print(f"Modello salvato con successo in: {model_output_path}")

    with open(args.out, 'w') as f:
        f.write('UniprotID,true_class,pred_class\n')
        for i in range(len(test_idx)):
            f.write(f'{df.iloc[test_idx[i]]["ID"]},{y_test.iloc[i]},{y_pred[i]}\n')

    if args.v=='y':
        with open(HERE/'SVM_performance.txt', 'w') as f:
            f.write(f'SVM \nHyperparameters: k={k}, C={c}, gamma={gamma}\n\n Model performance on benchmark set:\nMCC: {mcc:.3f}\nF1 Score: {f1:.3f}\nConfusion Matrix:\n{cm}\nPrecision: {precision:.3f}\nRecall: {recall:.3f}\nBalanced Accuracy: {bal_acc:.3f}\n\nFeatures ordered by RF importance: {features}\n')
        

