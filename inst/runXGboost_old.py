import sys
import os
import numpy as np 
import pandas as pd
import sklearn 
import xgboost as xgb
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold, train_test_split
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import matplotlib as mpl
from scipy import interp
from sklearn.preprocessing import scale
from sklearn.metrics import roc_auc_score, classification_report, accuracy_score, roc_curve, confusion_matrix, average_precision_score, precision_recall_curve
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold, train_test_split
from xgboost import XGBClassifier
import itertools
#import glmnet
import xgboost as xgb

from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

 

def runXGboost(labelData_path,predData_path,workingdir,dataName):
    os.chdir(workingdir)
    labelData = pd.read_csv(labelData_path,index_col=0)
    predData = pd.read_csv(predData_path,index_col=0)
    X,X_test, y, y_test = train_test_split(labelData.drop('label', axis=1), labelData['label'], test_size=0.33, random_state=42)
    #test_id=
    print(X)

    test=X_test
    test.to_csv('test.csv')
    y_test.to_csv('test_y.csv')
    i=1
    kf = StratifiedKFold(n_splits=5,random_state=1,shuffle=True)
    acc=[]
    for train_index,test_index in kf.split(X,y):
         print('\n{} of kfold {}'.format(i,kf.n_splits))
        # print(train_index)
         xtr,xvl = X.iloc[train_index],X.iloc[test_index]
         ytr,yvl = y.iloc[train_index],y.iloc[test_index]
         model= xgb.XGBClassifier(objective= 'binary:logistic')
         model.fit(xtr, ytr)
         pred=model.predict(xvl)
         print(pred)
         print(yvl)
         
         dfi=pd.DataFrame({'PRED': pred,'YVAL': yvl})
         print(dfi)
         fi='/Users/jyxi/Downloads/validation_'+dataName+str(i)+'.csv'
         dfi.to_csv(fi)
         #write.csv(dfi,paste0('/Users/jyxi/Downloads/validation_',i,'.csv')
         print('accuracy_score',accuracy_score(yvl,pred))
         acc=acc+[accuracy_score(yvl,pred)]
         i+=1
    print(sum(acc)/len(acc))




runXGboost(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])


