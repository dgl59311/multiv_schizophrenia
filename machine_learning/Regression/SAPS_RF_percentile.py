from sklearn.exceptions import ConvergenceWarning
import warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import pickle
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.base import clone
from time import time
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestRegressor
import copy
from scipy.stats import pearsonr

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

# details for RF
n_folds = 3
n_repeats = 20
random_state = 42
n_jobs = 9

model = Pipeline(steps=[('imputer', SimpleImputer(missing_values=-99999,
                                                  strategy='median')),
                        ('scale', RobustScaler()),
                        ('rf', RandomForestRegressor(random_state=random_state, n_estimators=1000, max_depth=10,
                                                     max_features='sqrt', n_jobs=n_jobs))])

# mean and standard deviation of the results in the test set
rmse_median_SAPS = []
rmse_25perc_SAPS = []
rmse_75perc_SAPS = []

r2_median_SAPS = []
r2_25perc_SAPS = []
r2_75perc_SAPS = []


r_median_SAPS = []
r_25perc_SAPS = []
r_75perc_SAPS = []


feature = []

# loop thru all the features
for count, ifeat in enumerate(Features.keys()):

    # append feature we are analysing
    feature.append(ifeat)

    # keep track of time it takes
    t1 = time()

    # load the feature we are interested in
    X = copy.deepcopy(Features[ifeat])
    X[np.isnan(X)] = -99999

    # for appending the rmse and r2 for each fold
    rmse_SAPS_test = []
    r2_SAPS_test = []
    r_SAPS_test = []

    # the split done here to make sure that all features have same splits
    cv_outer = RepeatedStratifiedKFold(n_splits=n_folds, n_repeats=n_repeats, random_state=random_state)

    # crossval for gridsearch
    cv_inner = n_folds

    # count outer split
    out_count = 0

    # loop through all folds and calculate the rmse and r2
    for train_indices, test_indices in cv_outer.split(X, y):

        # split into train data and test data
        X_train = X[train_indices, :]
        X_test = X[test_indices, :]

        SAPS_train = SAPS[train_indices]
        SAPS_test = SAPS[test_indices]

        # fit the data
        model_SAPS = clone(model)
        model_SAPS.fit(X_train, SAPS_train)

        # do prediction on the test data
        predicted_SAPS_test = model_SAPS.predict(X_test).ravel()
        predicted_SAPS_test[predicted_SAPS_test < 0] = 0
        predicted_SAPS_test[predicted_SAPS_test > 20] = 20

        # calculate metrics for the train data
        SAPS_rmse = np.sqrt(mean_squared_error(SAPS_test, predicted_SAPS_test))

        SAPS_r2 = r2_score(SAPS_test, predicted_SAPS_test)

        SAPS_r = pearsonr(SAPS_test, predicted_SAPS_test)[0]

        if np.isnan(SAPS_r):
            SAPS_r = 0

        # append the fold result
        rmse_SAPS_test.append(SAPS_rmse)
        r2_SAPS_test.append(SAPS_r2)
        r_SAPS_test.append(SAPS_r)

        # add count
        out_count += 1

    # append the mean and standard deviation across folds
    rmse_median_SAPS.append(np.median(rmse_SAPS_test))
    rmse_25perc_SAPS.append(np.percentile(rmse_SAPS_test, 25))
    rmse_75perc_SAPS.append(np.percentile(rmse_SAPS_test, 75))

    r2_median_SAPS.append(np.median(r2_SAPS_test))
    r2_25perc_SAPS.append(np.percentile(r2_SAPS_test, 25))
    r2_75perc_SAPS.append(np.percentile(r2_SAPS_test, 75))

    r_median_SAPS.append(np.median(r_SAPS_test))
    r_25perc_SAPS.append(np.percentile(r_SAPS_test, 25))
    r_75perc_SAPS.append(np.percentile(r_SAPS_test, 75))

    ## ----- print some stuff
    print("#---------Results--------------------------------")
    print(f"Total time feature {count + 1} '{ifeat}': {time() - t1}")
    print(f"SAPS median r2: {np.median(r2_SAPS_test)}")
    print(f"SAPS median r: {np.median(r_SAPS_test)}")
    print("#------------------------------------------------")

# create dictionary to save the results
results = {"feature" : feature,

           "rmse_median_SAPS" : rmse_median_SAPS,
           "rmse_25perc_SAPS" : rmse_25perc_SAPS,
           "rmse_75perc_SAPS" : rmse_75perc_SAPS,

           "r2_median_SAPS" : r2_median_SAPS,
           "r2_25perc_SAPS" : r2_25perc_SAPS,
           "r2_75perc_SAPS" : r2_75perc_SAPS,

           "r_median_SAPS" : r_median_SAPS,
           "r_25perc_SAPS" : r_25perc_SAPS,
           "r_75perc_SAPS" : r_75perc_SAPS}

# results in csv
results_df = pd.DataFrame.from_dict(results)

# save the results and be happy
results_df.to_csv("./20210810/SAPS_RF_percentile.csv")

# save the pickle files
with open("./20210810/SAPS_RF_percentile.pkl", "wb") as f:
    pickle.dump(results, f)