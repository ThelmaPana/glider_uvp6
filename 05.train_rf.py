#!/usr/bin/env python

'''
Train a RF on the training set.
Compute prediction metrics on the test set.
Predict new data, save predictions, which will be imported into Ecotaxa project 7587.
'''

import os
import glob
import pandas as pd
import math
import tarfile
import shutil
import datetime
import yaml
import numpy as np
from plotnine import *

import lib.read_settings as read_settings
import lib.datasets as datasets
import lib.rf as rf

# Read config
with open(r'config.yaml') as config_file:
    cfg = yaml.safe_load(config_file)
n_jobs = 24
random_state = 12

## Prepare training data
# Read learning set
df_learn = pd.read_csv('data/dataset/learn/04.learn_dataset.csv')

# Read test set
df_test = pd.read_csv('data/dataset/test/04.test_dataset.csv')
df_test_c = df_test.drop(df_test.loc[:, 'objid':'taxon'].columns, axis=1)

# Read new data
df_new = pd.read_csv('data/dataset/test/04.new_data.csv')
df_new_c = df_new.drop(df_new.loc[:, 'objid':'longitude'].columns, axis=1)

# Dataset composition
df_comp = df_learn.groupby(['taxon','set']).size().unstack(fill_value=0)
df_classes = df_learn[['taxon', 'plankton']].drop_duplicates().sort_values('taxon').reset_index(drop=True)
classes = np.array(df_classes['taxon'].tolist())

# Separate training from validation
df_train = df_learn[df_learn['set'] == 'train'].reset_index(drop = True).drop(df_learn.loc[:, 'objid':'set'].columns, axis=1)
df_valid = df_learn[df_learn['set'] == 'val'].reset_index(drop = True).drop(df_learn.loc[:, 'objid':'set'].columns, axis=1) 

# Generate class weights
class_counts = df_train.groupby('taxon').size()
count_max = 0
class_weights = {}
for idx in class_counts.items():
    count_max = (idx[1], count_max) [idx[1] < count_max]
for i,idx in enumerate(class_counts.items()):
    class_weights.update({idx[0] : math.sqrt(count_max / idx[1])})
    

## Gridsearch
gs_results, best_params = rf.gridsearch_rf(
    df_train, 
    df_valid, 
    classes = df_classes.taxon.tolist(),
    eval_metric = 'log_loss',
    max_features_try = [4, 6, 8, 10],
    min_samples_leaf_try = [2, 5, 10],
    n_estimators_try = [100, 200, 350, 500],
    output_dir = None,
    n_jobs = n_jobs,
    class_weights = class_weights,
    random_state = random_state
    )

# Plot gridsearch results
(ggplot(gs_results) +
  geom_point(aes(x='max_features', y='valid_log_loss', colour='factor(n_estimators)'))+
  facet_wrap('~min_samples_leaf', labeller = 'label_both') +
  labs(colour='n_estimators', title = 'Gridsearch results'))


## Train model with best hyperparameters
n_estimators = best_params['n_estimators']
max_features = best_params['max_features']
min_samples_leaf = best_params['min_samples_leaf']

# Fit the RF on training data
rf_fit = rf.train_rf(
    df = df_train, 
    n_estimators = n_estimators, 
    max_features = max_features, 
    min_samples_leaf = min_samples_leaf, 
    n_jobs = n_jobs, 
    class_weights = class_weights,
    random_state = random_state
)


## Predictions
# Test data
test_pred = rf_fit.predict_proba(df_test_c)

# Get classes and scores
df_test['y_pred'] = classes[np.argmax(test_pred, axis=1)]
df_test['score_pred'] = np.max(test_pred, axis=1)

# Save predictions
df_test.to_csv('data/dataset/test/05.test_dataset_predictions_rf.csv', index = False)


# New data
new_pred = rf_fit.predict_proba(df_new_c)

# Get classes and scores
df_new['y_pred'] = classes[np.argmax(new_pred, axis=1)]
df_new['score_pred'] = np.max(new_pred, axis=1)

# Save predictions
df_new.to_csv('data/dataset/test/05.new_data_predictions_rf.csv', index = False)


## Prepare predictions for import to Ecotaxa, but do not overwrite already validated objects
df_learn = pd.read_csv('data/dataset/learn/04.learn_dataset.csv')
#df_new = pd.read_csv('data/dataset/test/05.new_data_predictions_rf.csv')

# List already validated objects
val = df_learn['objid'].tolist()

# Extract predictions for non validated objects
df_new = df_new[~df_new['objid'].isin(val)].reset_index(drop = True)
df_new = df_new[['orig_id', 'y_pred', 'score_pred']]
df_new = df_new.rename(columns={
    'orig_id':'object_id',
    'y_pred':'object_annotation_category',
    'score_pred':'object_rf_score'
})


# Prepare formated tsv for Ecotaxa
df_new.loc[-1] = ['[t]', '[t]', '[f]']  # adding a row
df_new.index = df_new.index + 1  # shifting index
df_new = df_new.sort_index()  # sorting by index

# Save tsv
df_new.to_csv('data/dataset/learn/05.ecotaxa_rf_predictions.tsv', index=False, sep='\t', header=True)


