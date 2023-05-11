from sklearn.exceptions import ConvergenceWarning
import warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import pickle
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.base import clone
from time import time
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
import copy
from scipy.stats import pearsonr, median_abs_deviation

# read the data of all different features
with open("all_data_schizophrenia_194_eeg.pickle", "rb") as f:
    sz_data = pickle.load(f)

# check the keys
sz_data.keys()

# delete the things we don't want
SANS, SAPS = sz_data['SANS'], sz_data['SAPS']
Features = sz_data.copy()
del Features['SANS'], Features['SAPS']

# to keep proportions
y = pd.cut(SAPS.ravel(), bins=[-np.inf, 4., 6., 8., 10., 12., np.inf], labels=[1, 2, 3, 4, 5, 6])

# check if did what we wanted
Features.keys()

# details for elasticnet
n_folds = 3
n_repeats = 20
random_state = 42
n_jobs = 7

params = dict()
params['elastic__alpha'] = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
params['elastic__l1_ratio'] = np.arange(0.02, 1.02, 0.02)

model = Pipeline(steps=[('imputer', SimpleImputer(missing_values=-99999,
                                                  strategy='median')),
                        ('scale', RobustScaler()),
                        ('elastic', ElasticNet(random_state=random_state))])

# mean and standard deviation of the results in the test set
rmse_mean_SAPS = []
# rmse_mean_SANS = []
rmse_std_SAPS = []
# rmse_std_SANS = []
r2_mean_SAPS = []
# r2_mean_SANS = []
r2_std_SAPS = []
# r2_std_SANS = []
r_mean_SAPS = []
# r_mean_SANS = []
r_std_SAPS = []
# r_std_SANS = []
feature = []

# loop thru all the features
for count, ifeat in enumerate(Features.keys()):

    # if count == 3:
    #     break

    # append feature we are analysing
    feature.append(ifeat)

    # keep track of time it takes
    t1 = time()

    # load the feature we are interested in
    X = copy.deepcopy(Features[ifeat])
    X[np.isnan(X)] = -99999

    # for appending the rmse and r2 for each fold
    rmse_SAPS_test = []
    # rmse_SANS_test = []
    r2_SAPS_test = []
    # r2_SANS_test = []
    r_SAPS_test = []
    # r_SANS_test = []

    # the split done here to make sure that all features have same splits
    cv_outer = RepeatedStratifiedKFold(n_splits=n_folds, n_repeats=n_repeats, random_state=random_state)

    # crossval for gridsearch
    cv_inner = n_folds

    # grid search to find the optimal number of parameters
    grid = RandomizedSearchCV(model, params, scoring='neg_mean_squared_error', n_iter=100,
                              cv=cv_inner, refit=True, n_jobs=n_jobs, random_state=random_state)

    # count outer split
    out_count = 0

    # loop through all folds and calculate the rmse and r2
    for train_indices, test_indices in cv_outer.split(X, y):

        # split into train data and test data
        X_train = X[train_indices, :]
        X_test = X[test_indices, :]

        SAPS_train = SAPS[train_indices]
        SAPS_test = SAPS[test_indices]

        # SANS_train = SANS[train_indices]
        # SANS_test = SANS[test_indices]

        # fit the data
        model_SAPS = clone(grid)
        # model_SANS = clone(grid)
        model_SAPS.fit(X_train, SAPS_train)
        # model_SANS.fit(X_train, SANS_train)

        # do prediction on the test data
        predicted_SAPS_test = model_SAPS.predict(X_test).ravel()
        predicted_SAPS_test[predicted_SAPS_test < 0] = 0
        predicted_SAPS_test[predicted_SAPS_test > 20] = 20
        # predicted_SANS_test = model_SANS.predict(X_test).ravel()

        # calculate metrics for the train data
        SAPS_rmse = np.sqrt(mean_squared_error(SAPS_test, predicted_SAPS_test))
        # SANS_rmse = np.sqrt(mean_squared_error(SANS_test, predicted_SANS_test))

        SAPS_r2 = r2_score(SAPS_test, predicted_SAPS_test)
        # SANS_r2 = r2_score(SANS_test, predicted_SANS_test)

        SAPS_r = pearsonr(SAPS_test, predicted_SAPS_test)[0]
        # SANS_r = pearsonr(SANS_test, predicted_SANS_test)[0]

        if np.isnan(SAPS_r):
            SAPS_r = 0
        # if np.isnan(SANS_r):
        #     SANS_r = 0

        # append the fold result
        rmse_SAPS_test.append(SAPS_rmse)
        # rmse_SANS_test.append(SANS_rmse)
        r2_SAPS_test.append(SAPS_r2)
        # r2_SANS_test.append(SANS_r2)
        r_SAPS_test.append(SAPS_r)
        # r_SANS_test.append(SANS_r)

        # add count
        out_count += 1

    # append the mean and standard deviation across folds
    rmse_mean_SAPS.append(np.mean(rmse_SAPS_test))
    # rmse_mean_SANS.append(np.mean(rmse_SANS_test))
    rmse_std_SAPS.append(np.std(rmse_SAPS_test))
    # rmse_std_SANS.append(np.std(rmse_SANS_test))
    r2_mean_SAPS.append(np.mean(r2_SAPS_test))
    # r2_mean_SANS.append(np.mean(r2_SANS_test))
    r2_std_SAPS.append(np.std(r2_SAPS_test))
    # r2_std_SANS.append(np.std(r2_SANS_test))
    r_mean_SAPS.append(np.mean(r_SAPS_test))
    # r_mean_SANS.append(np.mean(r_SANS_test))
    r_std_SAPS.append(np.mean(r_SAPS_test))
    # r_std_SANS.append(np.std(r_SANS_test))

    ## ----- print some stuff
    print("#---------Results--------------------------------")
    print(f"Total time feature {count + 1} '{ifeat}': {time() - t1}")
    print(f"SAPS median r2: {np.mean(r2_SAPS_test)}")
    print(f"SAPS median r: {np.mean(r_SAPS_test)}")
    # print(f"SANS mean r2: {np.mean(r2_SANS_test)}")
    print("#------------------------------------------------")

# create dictionary to save the results
results = {"feature" : feature,
           "rmse_mean_SAPS" : rmse_mean_SAPS,
           # "rmse_mean_SANS" : rmse_mean_SANS,
           "rmse_std_SAPS" : rmse_std_SAPS,
           # "rmse_std_SANS" : rmse_std_SANS,
           "r2_mean_SAPS" : r2_mean_SAPS,
           # "r2_mean_SANS" : r2_mean_SANS,
           "r2_std_SAPS" : r2_std_SAPS,
           # "r2_std_SANS" : r2_std_SANS,
           "r_mean_SAPS" : r_mean_SAPS,
           # "r_mean_SANS" : r_mean_SANS,
           "r_std_SAPS" : r_std_SAPS}
           # "r_std_SANS" : r_std_SANS}

# results in csv
results_df = pd.DataFrame.from_dict(results)

# # save the results and be happy
# results_df.to_csv("ElasticNet_Nested_CrossValidation.csv")
#
# # save the pickle files
# with open("ElasticNet_Nested_CrossValidation.pkl", "wb") as f:
#     pickle.dump(results, f)