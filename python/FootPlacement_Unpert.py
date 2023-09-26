import os

# use utilities
import numpy as np
from scipy.interpolate import interp1d
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import pandas as pd
import csv
from utilities import *

# path information
# datapath = "C:/Users/mat950/Documents/Data/DataMoira_CSV/AnkleMoment"
datapath = "C:/Users/mat950/Documents/Data/DataMoira_CSV/FP_Force_COM"

# Create datatable for fit ?
BoolCreateTable = True

# init object to evaluate controller
myModel = RegressionModel(50000)

# time delay information
ths_delay = 0.2  # evaluate feedback model at 20% stance phase
tau = 0.1  # feedback delay in model

# selected subjects
#sVect = [3, 4, 7, 8, 14, 18, 21, 22, 24, 25, 31, 36, 41, 44]
sVect = [3]

# walking conditions
WalkCond = ["SteadyState_normal", "SteadyState_slow"]
# WalkCond = ['SteadyState_normal']

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

# settings
centerdata = True
removeorigin = True
pred_samples = range(0, 50)
order = 2
boolplot = True

# results output structure
Results = {}
for cond in WalkCond:
    Results[cond] = {}


if BoolCreateTable:
    # loop over all subjects
    for s in sVect:
        ctcond = 0
        print("started with subject" + str(s))
        for cond in WalkCond:
            ctcond = ctcond + 1
            # id of trials
            id = s * 10 + ctcond
            # check if datafile exists
            filepath_data = datapath + "/S" + str(s) + "/" + cond + "_data.csv"
            filepath_event = datapath + "/S" + str(s) + "/" + cond + "_event.csv"
            if os.path.exists(filepath_data) and os.path.exists(filepath_event):
                Dat = pd.read_csv(filepath_data)
                Event = pd.read_csv(filepath_event)

                # get inputs (mediolateral kinematics)
                t = np.array(Dat.time)
                Tankle = np.array(Dat.TAnkleLx)
                COM = np.array(Dat.COMx)
                FootL = np.array(Dat.FootLx)
                FootR = np.array(Dat.FootRx)

                # fitler input data
                #filtercutoff = 6
                #filteroder = 2
                #fs = 1 / np.nanmean(np.diff(t))
                #COM= ButterFilter_Low_NaNs(fs, COM, filteroder, filtercutoff)
                #FootL= ButterFilter_Low_NaNs(fs, FootL, filteroder, filtercutoff)
                #FootR= ButterFilter_Low_NaNs(fs, FootR, filteroder, filtercutoff)

                # run the foot placement model on this dataset
                output = FootPlacement_COMstate(
                    t,
                    COM,
                    FootL,
                    FootR,
                    Event,
                    centerdata,
                    removeorigin,
                    pred_samples,
                    order,
                    boolplot,
                    cond,
                )

                Results[cond][s] = output

# labels plot
ax1.set_xlabel("% single stance phase")
ax1.set_ylabel("Variance explained")
ax1.set_title(WalkCond[0])

ax2.set_xlabel("% single stance phase")
ax2.set_title(WalkCond[1])
plt.show()
