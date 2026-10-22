import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import matthews_corrcoef

RANDOM_STATE = 111


X=pd.read_csv('/home/um40/LB2_project_Group_1/Feature_Selection/feature.csv')

X.sort_values(by='id', inplace=True)

df=pd.read_csv('/home/um40/LB2_project_Group_1/Data_preparation/Cleaned_merged_dataset.tsv',sep='\t')

df = df.sort_values(by='ID')

y= df['class']
X= X.drop(columns='id', axis=1)
fold= df['fold']

feature_names=X.columns.to_list()




# K_GRID=[len(feature_names)//4,len(feature_names)//3,len(feature_names)//2,len(feature_names)]
K_GRID=[i for i in range(11,len(feature_names),2)]
C_GRID=[1,4,10,16]


GRID=[]
for k in K_GRID:
    for c in C_GRID:
        GRID.append((k,c,np.round(0.5/k,4))) 
        GRID.append((k,c,np.round(1/k,4))) 
        GRID.append((k,c,np.round(2/k,4))) 
        GRID.append((k,c,np.round(4/k,4))) 




def train_and_test(X_train, y_train, X_test, y_test, new_model, c, gamma):
    # train a classifier on (X_train, y_train) and return its mcc on (X_test, y_test)
    
    if X_train.shape[1] == 0:
        # no feature selected: the best we can do is predict the most frequent class
        majority = np.bincount(y_train).argmax()
        return np.mean(y_test == majority)
    
    # 1. learn the scaling (mean, std of each column) on the TRAINING data only
    scaler = StandardScaler().fit(X_train) # --> zeta scores = subtract the mean and divide by std (dependency from data, we compute mean & std on training deta => NOT USE TEST DATA TO COMPUTE THESE METRICS)
    # 2. train the classifier on the scaled training data
    model = new_model(c,gamma)
    model.fit(scaler.transform(X_train), y_train)
    # 3. scale the test data with the SAME mean and std, predict, measure mcc
    y_pred = model.predict(scaler.transform(X_test))
    return matthews_corrcoef(y_test, y_pred)



def choose_setting(select, grid, X, y, new_model,test_fold,fold):
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
             
            mcc = train_and_test(X.iloc[in_train, cols], y.iloc[in_train],
                                 X.iloc[in_val, cols],   y.iloc[in_val], new_model,c,gamma) #--> train a new model with the subset of best features, then test it on validation fold
            val_mccs.append(mcc)
        mean_val_mcc[setting] = np.mean(val_mccs)

    best = max(mean_val_mcc, key=mean_val_mcc.get)   # ties -> the first (smallest) one
    # if we have two combination of features with same mcc, always choose the one with fewer features (= less model's params, less overfit)
    return best, mean_val_mcc



results = []   # one entry per method

def nested_cv(name, family, select, grid, new_model,X,y,fold,verbose=True, store=True):
    
    start = time.time()
    test_mcc, chosen, selected_sets, inner_curves = [], [], [], []
    # selected_sets --> selected set of feature for each fold (?)

    for i in range(1,6):

        test_idx=fold[fold==i].index.to_list()
        train_idx=fold[(fold!=i) &( fold>0)].index.to_list()


        X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
        X_test,  y_test  = X.iloc[test_idx],  y.iloc[test_idx]

        # STEP 1 - inner loop: choose the setting using the training part only        

        best, curve = choose_setting(select, grid, X, y, new_model,i,fold) #--> choose best setting (= #best_features to keep) on TRAIN SET

        k, c, gamma = best

        # STEP 2 - redo the selection on the whole outer-training part
    
        cols = select(X_train, y_train, k) #--> given the #best_features (just a number!), select which are the #best most significant features on TRAIN TEST

        # STEP 3 - train on the outer-training part, test ONCE on the outer-test fold
        mcc = train_and_test(X_train.iloc[:, cols], y_train, X_test.iloc[:, cols], y_test, new_model, c, gamma) #--> train a new model with the subset of best features, then test it on validation fold

        test_mcc.append(mcc)
        chosen.append(best)
        selected_sets.append({feature_names[j] for j in cols})
        inner_curves.append(curve)
        if verbose:
            print(f"  outer fold {i}: chosen setting = {best}  "
                  f"features kept = {len(cols):2d}   test mcc = {mcc:.3f}")

    seconds = time.time() - start
    test_mcc = np.array(test_mcc)
    if verbose:
        print(f"{name}: {test_mcc.mean():.3f} ± {test_mcc.std():.3f}   ({seconds:.1f}s)\n")

    entry = dict(method=name, family=family, chosen=chosen,
                 k=[len(s) for s in selected_sets], mcc=test_mcc.mean(),
                 std=test_mcc.std(), seconds=seconds, fold_mccs=test_mcc,
                 sets=selected_sets, inner_curves=inner_curves, select=select, grid=grid)
    if store:
        results.append(entry)
    return entry



def select_rf(X_train, y_train, k):
    rf = RandomForestClassifier(n_estimators=100, random_state=0)
    rf.fit(X_train, y_train)
    ranking = np.argsort(rf.feature_importances_)[::-1]   # most important first (= gini importance)
    return ranking[:k]


def new_svm(c,gamma):
    return SVC(kernel="rbf", C=c, gamma=gamma)



if __name__=='__main__':
    svm_rf = nested_cv("Embedded: RF importance", "embedded", select_rf, GRID,new_svm,X,y,fold)
    with open('results_svm_rf.txt', 'w') as f:
        for k,v in svm_rf.items():
            print(f'{k}: {v}', file=f)

    