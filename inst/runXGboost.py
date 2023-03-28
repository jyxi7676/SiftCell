import sys
import os
import numpy as np 
import pandas as pd
import sklearn 
import xgboost as xgb
def runXGboost(labelData_path,predData_path,workingdir,dataName):
    os.chdir(workingdir)
    labelData = pd.read_csv(labelData_path,index_col=0)
    predData = pd.read_csv(predData_path,index_col=0)
    X_train = labelData.drop('label', axis=1)
    y_train = labelData['label']
    #Xgboost fitting
    clf = xgb.XGBClassifier(objective= 'binary:logistic')
    clf.fit(X_train,y_train)
    predData.columns= X_train.columns
    result = clf.predict(predData)
    unique, counts = np.unique(result, return_counts=True)
    positive_ind = labelData['label'] == "Nonsoup"
    positive=labelData[positive_ind]
    names1 = list(positive.index.values)
    positive_ind2 = result == "Nonsoup"
    names2=list(predData.index.values[positive_ind2])
    cellContaining = names1+names2 
    prefix = dataName+"Org_"
    cellContaining=[s.strip(prefix) for s in cellContaining]# remove the 8 from the string borders
    with open('CellContaining.txt', 'w') as f:
        for item in cellContaining:
            f.write("%s\n" % item)
    #print("Done!")
    #print (len(cellContaining)," Cell-containing droplets detected!")


runXGboost(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
