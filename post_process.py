
# coding: utf-8

# Call with either four arguments supplying the train and test input files and the locations of train and test output files
# or with seven arguments: four as above and then arguments specifying the sensitive attribute, score, and y if these are 
# different from attr, s, y


# In[1]:

import numpy as np
import pandas as pd
import sys
import cvxpy as cvx
from equalized_odds_and_calibration.eq_odds import Model


# In[2]:

train_file = sys.argv[1]
test_file = sys.argv[2]
train_out = sys.argv[3]
test_out = sys.argv[4]
if len(sys.argv) > 5:
    attr = sys.argv[5]
    score = sys.argv[6]
    y = sys.argv[7]
else: 
    attr = "attr"
    score = "s"
    y = "y"


# In[3]:

train = np.loadtxt(train_file, skiprows=1, delimiter = ",")
test = np.loadtxt(test_file, skiprows=1, delimiter = ",")


# In[4]:

cols = pd.read_csv(train_file).columns
train_A = train[train[:, cols.get_loc(attr)] ==1]
train_B = train[train[:, cols.get_loc(attr)] ==0]
test_A = test[test[:, cols.get_loc(attr)] ==1]
test_B = test[test[:, cols.get_loc(attr)] ==0]


# # Hardt post-processing

# # observed

# In[5]:

model_A = Model(train_A[:, cols.get_loc(score)], train_A[:, cols.get_loc(y)])
model_B = Model(train_B[:, cols.get_loc(score)], train_B[:, cols.get_loc(y)])
model_test_A = Model(test_A[:, cols.get_loc(score)], test_A[:, cols.get_loc(y)])
model_test_B = Model(test_B[:, cols.get_loc(score)], test_B[:, cols.get_loc(y)])


# In[6]:

# Find mixing rates for post-processing method
_, _, mix_rates = Model.eq_odds(model_A, model_B)
model_corrected_A, model_corrected_B = Model.eq_odds(model_A, model_B, mix_rates)
model_corrected_test_A, model_corrected_test_B = Model.eq_odds(model_test_A, model_test_B, mix_rates)


# In[7]:

test_B = np.append(test_B, np.expand_dims(model_corrected_test_B.pred, axis =1), 1)
test_A = np.append(test_A, np.expand_dims(model_corrected_test_A.pred, axis =1), 1)
train_B = np.append(train_B, np.expand_dims(model_corrected_B.pred, axis =1), 1)
train_A = np.append(train_A, np.expand_dims(model_corrected_A.pred, axis =1), 1)


# In[8]:

cols = np.concatenate((np.array(cols), np.array(['eo_fair_pred'])))


# In[9]:

test_B_pd = pd.DataFrame(data = test_B, columns = cols)
test_A_pd = pd.DataFrame(data = test_A, columns = cols)
train_B_pd = pd.DataFrame(data = train_B, columns = cols)
train_A_pd = pd.DataFrame(data = train_A, columns = cols)
train = pd.concat([train_A_pd, train_B_pd], axis = 0)
test = pd.concat([test_A_pd, test_B_pd], axis = 0)


# In[10]:

test.to_csv(test_out)
train.to_csv(train_out)




