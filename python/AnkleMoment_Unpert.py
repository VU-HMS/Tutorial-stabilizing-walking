import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


# use utilities
from utilities import *

# path information
datapath = "C:/Users/mat950/Documents/Data/DataMoira_CSV/AnkleMoment"

# Create datatable for fit ?
BoolCreateTable = True

# init object to evaluate controller
myModel = RegressionModel(100000)

# time delay information
ths_delay = 0.3  # evaluate feedback model at 20% stance phase
tau = 0.1  # feedback delay in model
ths_delay_vect = np.array(
    [
        0.05,
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
    ]
)

# selected subjects
sVect = [3, 7, 8, 14, 18, 22, 24, 25, 31, 36, 41, 44]

# walking conditions
WalkCond = ["SteadyState_normal", "SteadyState_slow"]

if BoolCreateTable:
    # loop over all subjects
    for s in sVect:
        ctcond = 0
        print("started with subject" + str(s))
        for cond in WalkCond:
            ctcond = ctcond + 1
            # id of trials
            # id = s*10 + ctcond
            # check if datafile exists
            filepath_data = datapath + "/S" + str(s) + "/" + cond + "_data.csv"
            filepath_event = datapath + "/S" + str(s) + "/" + cond + "_event.csv"
            if os.path.exists(filepath_data) and os.path.exists(filepath_event):
                Dat = pd.read_csv(filepath_data)
                Event = pd.read_csv(filepath_event)
                # extract relevant data from .csv file
                t = np.array(Dat.time)
                Tankle = np.array(Dat.TAnkleLy)
                COM = np.array(Dat.COMy)
                FootL = np.array(Dat.FootLy)
                # filter input data
                filtercutoff = 6
                filteroder = 2
                fs = 1 / np.nanmean(np.diff(t))
                COMfilt = ButterFilter_Low_NaNs(fs, COM, filteroder, filtercutoff)
                COMdfilt = central_difference(t, COMfilt)
                Tanklefilt = ButterFilter_Low_NaNs(
                    fs, Tankle, filteroder, filtercutoff
                )
                FootLfilt = ButterFilter_Low_NaNs(
                    fs, FootL, filteroder, filtercutoff
                )

                # nondim all inputs
                g = 9.81
                m = np.mean(Dat.GRFLz + Dat.GRFRz) / g
                L = np.mean(Dat.COMz)
                Tankle = Tanklefilt / (m * L * g)
                r_COM_Foot = (COMfilt - FootLfilt) / L
                COMd = COMdfilt / np.sqrt(g * L)
                # loop over left heelstrikes
                n_lhs = len(Event.lhs)
                ctDelay = -1
                for iDelay in ths_delay_vect:
                    ctDelay = ctDelay + 1
                    id = s * 10 + ctcond + ctDelay * 1000
                    ths_delay = iDelay
                    for i in range(1, n_lhs - 3):
                        # find toe-off
                        iabove = Event.lto > Event.lhs[i]
                        dt_stance = Event.lto[np.argmax(iabove)] - Event.lhs[i]
                        # get COM state information
                        iabove = t > (Event.lhs[i] + ths_delay * dt_stance)
                        ix = np.argmax(iabove)
                        COM_input = COM[ix]
                        COMd_input = COMd[ix]
                        # get Torque information
                        iabove = t > (Event.lhs[i] + ths_delay * dt_stance + tau)
                        iy = np.argmax(iabove)
                        Tankle_input = Tankle[iy]
                        # add datapoint
                        if (
                            not (np.isnan(COM_input))
                            and not (np.isnan(COMd_input))
                            and not (np.isnan(Tankle_input))
                        ):
                            myModel.adddatapoint(
                                COM_input, COMd_input, Tankle_input, id
                            )

    myModel.writetocsv("RegressionModel_DataMoira.csv")
else:
    # init object to evaluate controller
    myModel = RegressionModel(50000)
    myModel.importfromcsv("RegressionModel_DataMoira.csv")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# plot figure of slow walking
RsqAll_slow = np.zeros(shape=(len(sVect), len(ths_delay_vect)))
ctDelay = -1
ctcond = 1
cti = 0
for i in ths_delay_vect:
    idVect = []
    ctDelay = ctDelay + 1
    cts = 0
    for s in sVect:
        id = s * 10 + ctcond + ctDelay * 1000
        idVect.append(id)
        myModel.fitmodel(id, False)  # fit model on individual
        RsqAll_slow[cts, cti] = myModel.r2
        cts = cts + 1
    cti = cti + 1

cts = 0
for s in sVect:
    # select datapoints of this subject
    ax1.plot(ths_delay_vect, RsqAll_slow[cts, :], color=(0.6, 0.6, 0.6), linewidth=1)
    cts = cts + 1
RsqMean = np.mean(RsqAll_slow, axis=0)
ax1.plot(ths_delay_vect, RsqMean, color=(0.2, 0.2, 0.2), linewidth=2)

# plot figure of normal walking speed => we need this one
RsqAll_fast = np.zeros(shape=(len(sVect), len(ths_delay_vect)))
ctDelay = -1
ctcond = 2
cti = 0
for i in ths_delay_vect:
    idVect = []
    ctDelay = ctDelay + 1
    cts = 0
    for s in sVect:
        id = s * 10 + ctcond + ctDelay * 1000
        idVect.append(id)
        myModel.fitmodel(id, False)  # fit model on individual
        RsqAll_fast[cts, cti] = myModel.r2
        cts = cts + 1
    cti = cti + 1

cts = 0
for s in sVect:
    # select datapoints of this subject
    ax2.plot(ths_delay_vect, RsqAll_fast[cts, :], color=(0.6, 0.6, 0.6), linewidth=1)
    cts = cts + 1
RsqMean = np.mean(RsqAll_fast, axis=0)
ax2.plot(ths_delay_vect, RsqMean, color=(0.2, 0.2, 0.2), linewidth=2)

ax1.set_ylim([0, 0.8])
ax2.set_ylim([0, 0.8])


# new figure for the tutorial paper
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
cts = 0
for s in sVect:
    # select datapoints of this subject
    ax.plot(ths_delay_vect, RsqAll_fast[cts, :], color=(0.6, 0.6, 0.6), linewidth=1)
    cts = cts + 1
RsqMean = np.mean(RsqAll_fast, axis=0)
ax.plot(ths_delay_vect, RsqMean, color=(0.2, 0.2, 0.2), linewidth=2)

for s in sVect:
    ctcond = 2
    ctDelay = 4
    id = s * 10 + ctcond + ctDelay * 1000
    idVect.append(id)
    myModel.fitmodel(id, True)  # fit model on individual

# show all figures
plt.show()

# print('Rsq fast walking '+ str(np.mean(RsqVect_fast)) +'  std ' + str(np.std(RsqVect_fast)))
# print('Rsq fast walking pooled '+ str(myModel.r2))
# plt.show()
