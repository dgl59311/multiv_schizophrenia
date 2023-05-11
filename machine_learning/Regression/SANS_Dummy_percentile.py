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
from sklearn.dummy import DummyRegressor
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
y = pd.cut(SANS.ravel(), bins=[-np.inf, 4., 6., 8., 10., 13., 16., np.inf], labels=[1, 2, 3, 4, 5, 6, 7])

# check if did what we wanted
Features = ["Dummy"]

# details for RF
n_folds = 3
n_repeats = 20
random_state = 42
n_jobs = -1

model = Pipeline(steps=[('imputer', SimpleImputer(missing_values=-99999,
                                                  strategy='median')),
                        ('regressor', DummyRegressor(strategy='median'))])

# mean and standard deviation of the results in the test set
rmse_median_SANS = []
rmse_25perc_SANS = []
rmse_75perc_SANS = []

r2_median_SANS = []
r2_25perc_SANS = []
r2_75perc_SANS = []


r_median_SANS = []
r_25perc_SANS = []
r_75perc_SANS = []


feature = []

# loop thru all the features
for count, ifeat in enumerate(Features):

    # append feature we are analysing
    feature.append(ifeat)

    # keep track of time it takes
    t1 = time()

    # load the feature we are interested in
    X = copy.deepcopy(sz_data["Activity"])
    X[np.isnan(X)] = -99999

    # for appending the rmse and r2 for each fold
    rmse_SANS_test = []
    r2_SANS_test = []
    r_SANS_test = []

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

        SANS_train = SANS[train_indices]
        SANS_test = SANS[test_indices]

        # fit the data
        model_SANS = clone(model)
        model_SANS.fit(X_train, SANS_train)

        # do prediction on the test data
        predicted_SANS_test = model_SANS.predict(X_test).ravel()
        predicted_SANS_test[predicted_SANS_test < 0] = 0
        predicted_SANS_test[predicted_SANS_test > 25] = 25

        # calculate metrics for the train data
        SANS_rmse = np.sqrt(mean_squared_error(SANS_test, predicted_SANS_test))

        SANS_r2 = r2_score(SANS_test, predicted_SANS_test)

        SANS_r = pearsonr(SANS_test, predicted_SANS_test)[0]

        if np.isnan(SANS_r):
            SANS_r = 0

        # append the fold result
        rmse_SANS_test.append(SANS_rmse)
        r2_SANS_test.append(SANS_r2)
        r_SANS_test.append(SANS_r)

        # add count
        out_count += 1

    # append the mean and standard deviation across folds
    rmse_median_SANS.append(np.median(rmse_SANS_test))
    rmse_25perc_SANS.append(np.percentile(rmse_SANS_test, 25))
    rmse_75perc_SANS.append(np.percentile(rmse_SANS_test, 75))

    r2_median_SANS.append(np.median(r2_SANS_test))
    r2_25perc_SANS.append(np.percentile(r2_SANS_test, 25))
    r2_75perc_SANS.append(np.percentile(r2_SANS_test, 75))

    r_median_SANS.append(np.median(r_SANS_test))
    r_25perc_SANS.append(np.percentile(r_SANS_test, 25))
    r_75perc_SANS.append(np.percentile(r_SANS_test, 75))

    ## ----- print some stuff
    print("#---------Results--------------------------------")
    print(f"Total time feature {count + 1} '{ifeat}': {time() - t1}")
    print(f"SANS median r2: {np.median(r2_SANS_test)}")
    print(f"SANS median r: {np.median(r_SANS_test)}")
    print("#------------------------------------------------")

# create dictionary to save the results
results = {"feature" : feature,

           "rmse_median_SANS" : rmse_median_SANS,
           "rmse_25perc_SANS" : rmse_25perc_SANS,
           "rmse_75perc_SANS" : rmse_75perc_SANS,

           "r2_median_SANS" : r2_median_SANS,
           "r2_25perc_SANS" : r2_25perc_SANS,
           "r2_75perc_SANS" : r2_75perc_SANS,

           "r_median_SANS" : r_median_SANS,
           "r_25perc_SANS" : r_25perc_SANS,
           "r_75perc_SANS" : r_75perc_SANS}

# results in csv
results_df = pd.DataFrame.from_dict(results)

# save the results and be happy
results_df.to_csv("./20210810/SANS_dummy_percentile.csv")

# save the pickle files
with open("./20210810/SANS_dummy_percentile.pkl", "wb") as f:
    pickle.dump(results, f)