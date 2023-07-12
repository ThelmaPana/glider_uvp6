import pandas as pd
import numpy as np
import itertools as it
import os
import pickle
import glob

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, log_loss

def gridsearch_rf(df_train, df_valid, classes, eval_metric, max_features_try, min_samples_leaf_try, n_estimators_try, output_dir, n_jobs, class_weights=None, random_state=None):
    """
    Perform a grid search to find best hyperparameters for random forest model.
    
    Args:
        df_train (DataFrame): training data to use to fit grid search
        df_valid (DataFrame): validation data to use to evaluate grid search
        classes (list, array): name of classes
        eval_metric (str): metric to use for hyperparameters selection ('accuracy', 'balanced_accuracy' or 'log_loss')
        max_features_try (list): tries for number of variables per node; default sqrt(nb of vars)
        min_samples_leaf_try (list): tries for min number of objects in leaf; default for classif = 5
        n_estimators_try (list): tries for number of estimators (usually between 100 and 500)
        output_dir (str): directory where to save gridsearch results
        n_jobs (int): number of cores to use 
        class_weights (dict): weights for classes
        random_state (int or RandomState): controls both the randomness of the bootstrapping and features sampling; default=None
    
    Returns:
        results (DataFrame): results of grid search
        best_params (dict): best parameters based on evaluation metric
    """
    
    # Shuffle data
    df_train = df_train.sample(frac=1, random_state=random_state).reset_index(drop=True)
    df_valid = df_valid.sample(frac=1, random_state=random_state).reset_index(drop=True)
    
    # Split data and labels
    y_train = df_train['taxon']
    X_train = df_train.drop('taxon', axis=1)
    
    y_valid = df_valid['taxon']
    X_valid = df_valid.drop('taxon', axis=1)
    
    # Prepare one-hot encoding to compute log-loss for validation data
    mlb = MultiLabelBinarizer(classes=classes)
    y_true = mlb.fit_transform([[l] for l in y_valid])
    
    # Build grid of hyperparameters to explore
    grid = {
        'max_features': max_features_try, 
        'min_samples_leaf': min_samples_leaf_try,
    }
    # Make a list of all parameters combinations
    keys = list(grid.keys())
    grid_list = list(it.product(*(grid[key] for key in keys)))

    # Initiate empty dict for results
    results = {
        'n_estimators': [],
        'max_features': [],
        'min_samples_leaf': [],
        'valid_accuracy': [],
        'valid_balanced_accuracy':[],
        'valid_log_loss':[]
    }

    # First loop on parameters other than n_estimators
    for max_features, min_samples_leaf in grid_list:
        print(f"Trying parameters max_features = {max_features} and min_samples_leaf = {min_samples_leaf}.")
    
        # Initiate a RF model with warm start
        rf = RandomForestClassifier(
            criterion='entropy', 
            min_samples_split=2, 
            max_features=max_features, 
            min_samples_leaf=min_samples_leaf,
            warm_start=True,
            n_jobs=n_jobs,
            class_weight=class_weights,
            random_state=random_state
        )
        
        # Second loop on n_estimators
        for n_estimators in n_estimators_try:
            print(f"Number of estimators = {n_estimators}")
            
            # Set number of estimators in RF model
            rf.n_estimators = n_estimators
            
            # Fit on training data
            rf.fit(X=X_train, y=y_train)
            
            # Compute accuracy on validation data
            valid_accuracy = accuracy_score(y_valid, rf.predict(X_valid))
            valid_balanced_accuracy = balanced_accuracy_score(y_valid, rf.predict(X_valid))
            
            # Compute log loss on validation data  
            # log_loss only accepts weights as sample_weights and not as class_weights, compute sample_weights
            if class_weights is not None:
                sample_weight = [class_weights[c] for c in y_valid]
            else:
                sample_weight = None
            valid_log_loss = log_loss(y_true, rf.predict_proba(X_valid), sample_weight=sample_weight)

            # Store results in dict
            results['n_estimators'].append(n_estimators)
            results['max_features'].append(max_features)
            results['min_samples_leaf'].append(min_samples_leaf)
            results['valid_accuracy'].append(valid_accuracy) 
            results['valid_balanced_accuracy'].append(valid_balanced_accuracy) 
            results['valid_log_loss'].append(valid_log_loss) 

    # Write gridsearch results
    #with open(os.path.join(output_dir, 'train_results.pickle'),'wb') as results_file:
    #    pickle.dump(results, results_file)
        
    # Convert to datfarame
    results = pd.DataFrame(results)
    
    # Extract best parameters based on evaluation metric value on validation data
    if eval_metric == 'log_loss':
        # if evaluation metric is log loss, look for the smallest value
        best_params = results.nsmallest(1, 'valid_'+ eval_metric).reset_index(drop=True).drop(['valid_accuracy', 'valid_balanced_accuracy', 'valid_log_loss'], axis=1)
    else:
        # in other cases, look for the largest value
        best_params = results.nlargest(1, 'valid_'+ eval_metric).reset_index(drop=True).drop(['valid_accuracy', 'valid_balanced_accuracy', 'valid_log_loss'], axis=1)
    best_params = best_params.iloc[0].to_dict()
        
    # Print GS results
    print(f"Hyperparameters selection using {eval_metric} value.")
    print(f"Selected hyperparameters are: n_estimators = {best_params['n_estimators']}, max_features = {best_params['max_features']}, min_samples_leaf = {best_params['min_samples_leaf']}")

    return results, best_params


def train_rf(df, n_estimators, max_features, min_samples_leaf, n_jobs, class_weights, random_state=None):
    """
    Fit a random forest model on data.
    
    Args:
        df (DataFrame): data to use for training
        tree_nb (int): number of trees for the RF model
        max_features (int): number of variables per node
        min_samples_leaf (int): min number of objects in leaf
        n_jobs (int): number of cores to use 
        class_weights(dict): weights for classes
        random_state (int or RandomState): controls both the randomness of the bootstrapping and features sampling; default=None
    
    Returns:
        rf (RandomForestClassifier): fitted random forest model
    """
    
    # Shuffle data
    df = df.sample(frac=1).reset_index(drop=True)
    
    # Split data and labels
    y_train = df['taxon']
    X_train = df.drop('taxon', axis=1)
    
    # Initiate RF model
    rf = RandomForestClassifier(
        n_estimators=n_estimators, 
        criterion='entropy', 
        min_samples_split=2, 
        min_samples_leaf=min_samples_leaf, 
        max_features=max_features,
        n_jobs=n_jobs,
        class_weight=class_weights,
        random_state=random_state
    )
    
    # Fit the RF model
    rf = rf.fit(X=X_train, y=y_train)
    
    return rf


def predict_evaluate_rf(rf_model, df, df_classes, output_dir):
    """
    Evaluate a random forest model.
    
    Args:
        rf_model (RandomForestClassifier): random forest model to evaluate
        df (DataFrame): data to use for model evaluation
        df_classes (DataFrame): dataframe of classes with living attribute
        output_dir (str): directory where to save prediction results
    
    Returns:
        nothing
    """
    
    # Split data and labels
    y = df['taxon'].tolist()
    y = np.array(y)
    X = df.drop('taxon', axis=1)

    # Make a list of classes
    classes = df_classes['classif_id'].tolist()
    classes.sort()
    classes = np.array(classes)
    
    # and of regrouped classes
    classes_g = df_classes['classif_id_2'].tolist()
    classes_g = list(set(classes_g))
    classes_g.sort()
    classes_g = np.array(classes_g)
    
    # Make a list of plankton classes
    plankton_classes = df_classes[df_classes['plankton']]['classif_id'].tolist()
    plankton_classes = np.array(plankton_classes)
    
    # Make a list of plankton classes for grouped classes
    plankton_classes_g = df_classes[df_classes['plankton']]['classif_id_2'].tolist()
    plankton_classes_g = np.array(plankton_classes_g)
    
    # Predict test data
    y_pred = rf_model.predict(X)
    
    # Compute accuracy between true labels and predicted labels
    accuracy = accuracy_score(y, y_pred)
    balanced_accuracy = balanced_accuracy_score(y, y_pred)
    plankton_precision = precision_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    plankton_recall = recall_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    
    # Display results
    print(f'Test accuracy = {accuracy}')
    print(f'Balanced test accuracy = {balanced_accuracy}')
    print(f'Weighted plankton precision = {plankton_precision}')
    print(f'Weighted plankton recall = {plankton_recall}')
     
    ## Now do the same after regrouping objects to larger classes
    # Generate taxonomy match between taxo used for classif and larger ecological classes 
    taxo_match = df_classes.set_index('classif_id').to_dict('index')
    
    # Convert true classes to larger ecological classes
    y_g = np.array([taxo_match[t]['classif_id_2'] for t in y])
    
    # Convert predicted classes to larger ecological classes
    y_pred_g = np.array([taxo_match[p]['classif_id_2'] for p in y_pred])
    
    # Compute accuracy, precision and recall for living classes and loss from true labels and predicted labels
    accuracy_g = accuracy_score(y_g, y_pred_g)
    balanced_accuracy_g = balanced_accuracy_score(y_g, y_pred_g)
    plankton_precision_g = precision_score(y_g, y_pred_g, labels=plankton_classes_g, average='weighted', zero_division=0)
    plankton_recall_g = recall_score(y_g, y_pred_g, labels=plankton_classes_g, average='weighted', zero_division=0)
    
    # Display results
    print(f'Grouped test accuracy = {accuracy_g}')
    print(f'Grouped balanced test accuracy = {balanced_accuracy_g}')
    print(f'Grouped weighted plankton precision = {plankton_precision_g}')
    print(f'Grouped weighted plankton recall = {plankton_recall_g}')

    # Write classes and test metrics into a test file
    with open(os.path.join(output_dir, 'test_results.pickle'),'wb') as test_file:
        pickle.dump({
            'classes': classes,
            'classes_g': classes_g,
            'plankton_classes': plankton_classes,
            'plankton_classes_g': plankton_classes_g,
            'true_classes': y,
            'predicted_classes': y_pred,
            'true_classes_g': y_g,
            'predicted_classes_g': y_pred_g,
            'accuracy': accuracy,
            'balanced_accuracy': balanced_accuracy,
            'plankton_precision': plankton_precision,
            'plankton_recall': plankton_recall,
            'accuracy_g': accuracy_g,
            'balanced_accuracy_g': balanced_accuracy_g,
            'plankton_precision_g': plankton_precision_g,
            'plankton_recall_g': plankton_recall_g,
        },
        test_file)

    pass


def predict_rf(rf_model, df, df_classes):
    """
    Evaluate a random forest model.
    
    Args:
        rf_model (RandomForestClassifier): random forest model to evaluate
        df (DataFrame): data to use for model evaluation
        df_classes (DataFrame): dataframe of classes with living attribute
        output_dir (str): directory where to save prediction results
    
    Returns:
        nothing
    """
    
    # Split data and labels
    y = df['taxon'].tolist()
    y = np.array(y)
    X = df.drop('taxon', axis=1)

    # Make a list of classes
    classes = df_classes['taxon'].tolist()
    classes.sort()
    classes = np.array(classes)
    
    # Make a list of plankton classes
    plankton_classes = df_classes[df_classes['plankton']]['taxon'].tolist()
    plankton_classes = np.array(plankton_classes)
    
    # Predict test data
    y_pred = rf_model.predict(X)
    
    # Compute accuracy between true labels and predicted labels
    accuracy = accuracy_score(y, y_pred)
    balanced_accuracy = balanced_accuracy_score(y, y_pred)
    plankton_precision = precision_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    plankton_recall = recall_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    
    # Display results
    print(f'Test accuracy = {accuracy}')
    print(f'Balanced test accuracy = {balanced_accuracy}')
    print(f'Weighted plankton precision = {plankton_precision}')
    print(f'Weighted plankton recall = {plankton_recall}')
     
    # Return predictions
    return y_pred
