# -*- coding: utf-8 -*-
import os 
import numpy as np
import pandas as pd
import matplotlib 
import sklearn as sk
from sklearn import cross_validation
from matplotlib import pyplot as plt
from sklearn import *
from funk_py import *
from sklearn.ensemble.forest import RandomForestRegressor

## Fonction pour que python reconnaisse les floats à virgules.
def funk_to_float(col):
    # col = array of dataframe column
    return([float(x.replace(",",".")) for x in col])

## Fonction pour spérarer les char des nums
def funk_data_cat_char_from_num(data):
    # data : data.frame
    col_tr_type = [type(data.loc[1,x]) for x in data.columns]
    col_tr_str = data.columns[[x == str for x in col_tr_type]]
    for esc in col_tr_str:
        try:
            data.loc[:,esc] = funk_to_float(data.loc[:,esc])
        except:
            print esc
            #print "\n"
            #print data.loc[:,esc].head()
    col_tr_type = [type(data.loc[1,x]) for x in data.columns]
    col_tr_str = data.columns[[x == str for x in col_tr_type]]
    test  = set(data.columns)
    return([data.loc[:,col_tr_str],data.loc[:,test.difference(set(col_tr_str))]])

#Fonctions de plits
def funk_split(data,ratio):
    # data : dataframe, ratio : float
    data_train, data_test = cross_validation.train_test_split(data, test_size=1-ratio, random_state=0)
    data_train = pd.DataFrame(data_train,columns=data.columns)
    data_test = pd.DataFrame(data_test,columns=data.columns)
    return([data_train,data_test])

# Fonction d'incrément pour les char
def funk_increment_col(data):
    # data : dataframe
    data_res = data.copy()
    for esc in data_res.columns:
        unik = data_res.loc[:,esc].unique()
        nunik = len(data_res.loc[:,esc].unique())
        for it in range(0,nunik):
            data_res.loc[data_res.loc[:,esc] == unik[it],esc] = it
    return(data_res)

# Fonction rechargeant les données spliter:
def funk_charge_set_data(strin_adress):
    # strin_adress : string
    l=list()    
    l.append(pd.read_csv(strin_adress+"\\data_char_train.csv"))
    l.append(pd.read_csv(strin_adress+"\\data_char_test.csv"))
    l.append(pd.read_csv(strin_adress+"\\data_num_train.csv"))
    l.append(pd.read_csv(strin_adress+"\\data_num_test.csv"))
    l.append(pd.read_csv(strin_adress+"\\trgt1_train.csv"))
    l.append(pd.read_csv(strin_adress+"\\trgt1_test.csv"))
    l.append(pd.read_csv(strin_adress+"\\trgt2_train.csv"))
    l.append(pd.read_csv(strin_adress+"\\trgt2_test.csv"))
    return(l)
# Fonction calculant l'entropie de Shanon
def funk_enthropy(test):
    # test : matrix or pandas.Dataframe  
    return(-sum([x*np.log(x) for x in test if x!=0]))













