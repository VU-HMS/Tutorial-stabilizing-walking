import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.signal import butter, filtfilt

# use utilities
from utilities import *

# path information
# datapath = "C:/Users/mat950/Documents/Data/DataMoira_CSV/AnkleMoment"
datapath = "C:/Users/mat950/Documents/Data/DataMoira_CSV/FP_Force_COM"
# Create datatable for fit ?
BoolCreateTable = True

# selected subjects
sVect = [3, 4, 7, 8, 14, 18, 21, 22, 24, 25, 36, 41, 44]
#sVect = [31]

# walking conditions
WalkCond = ["SteadyState_normal", "SteadyState_slow"]

# pre allocate ouput vars
nsubj = len(sVect)
nConds = len(WalkCond)
corr_phase_xcom = np.full([50, 2, 2, nConds, nsubj], np.nan)
gain_phase_xcom = np.full([50, 2, 2, nConds, nsubj], np.nan)
lag_xcom = np.full([2, 2, nConds, nsubj], np.nan)
corr_phase_com = np.full([50, 2, 2, nConds, nsubj], np.nan)
gain_phase_com = np.full([50, 2, 2, nConds, nsubj], np.nan)
gain_phase_vcom = np.full([50, 2, 2, nConds, nsubj], np.nan)
lag_com = np.full([2, 2, nConds, nsubj], np.nan)

if BoolCreateTable:
    # loop over all subjects
    cts = -1
    for s in sVect:
        cts = cts + 1
        ctcond = -1
        print("started with subject" + str(s))
        for cond in WalkCond:
            ctcond = ctcond + 1
            # check if datafile exists
            filepath_data = datapath + "/S" + str(s) + "/" + cond + "_data.csv"
            filepath_event = datapath + "/S" + str(s) + "/" + cond + "_event.csv"
            if os.path.exists(filepath_data) and os.path.exists(filepath_event):
                Dat = pd.read_csv(filepath_data)
                Event = pd.read_csv(filepath_event)
                # extract input information
                COM = np.vstack((Dat.COMx, Dat.COMy)).T
                time = np.array(Dat.time).T

                FootL = np.vstack((Dat.FootLx, Dat.FootLy)).T
                FootR = np.vstack((Dat.FootRx, Dat.FootRy)).T
                GRFR = np.vstack((Dat.GRFRx, Dat.GRFRy)).T
                GRFL = np.vstack((Dat.GRFLx, Dat.GRFLy)).T

                # filter COM data
                fs = 1 / np.nanmean(np.diff(time))
                COM_filt = ButterFilter_Low_NaNs(fs, COM, 2, 6)
                COMd_filt = central_difference(time, COM_filt)

                # filter GRF dat
                GRFR_filt = ButterFilter_Low_NaNs(fs, GRFR, 2, 6)
                GRFL_filt = ButterFilter_Low_NaNs(fs, GRFL, 2, 6)
                GRF = GRFR_filt + GRFL_filt

                # filter foot positions
                FootL_filt = ButterFilter_Low_NaNs(fs, FootL, 2, 6)
                FootR_filt = ButterFilter_Low_NaNs(fs, FootR, 2, 6)

                # compute extrapolated center of mass
                L = np.nanmean(Dat.COMz)
                xCOM = COM_filt + COMd_filt / np.sqrt(9.81 / L)

                # relate xCOM information and GRF information
                t_lhs = np.array(Event.lhs)
                t_rhs = np.array(Event.rhs)
                maxlag = 40
                R = feedback_com_xcom(
                    GRF, COM_filt, COMd_filt, xCOM, FootL_filt, FootR_filt, 40, t_lhs, t_rhs, time
                )

                # store results in nd arrays
                corr_phase_xcom[:, :, :, ctcond, cts] = R[0]
                gain_phase_xcom[:, :, :, ctcond, cts] = R[1]
                lag_xcom[:, :, ctcond, cts] = R[2]
                corr_phase_com[:, :, :, ctcond, cts] = R[3]
                gain_phase_com[:, :, :, ctcond, cts] = R[4]
                gain_phase_vcom[:, :, :, ctcond, cts] = R[5]
                lag_com[:, :, ctcond, cts] = R[6]


print("To Do: investigate why results are slightly different")
fig = plt.figure()

x = np.linspace(51, 100, 50)
ax1 = fig.add_subplot(2, 2, 1)
dim = 0  # x-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = corr_phase_xcom[:, dim, ft, ctcond, :]
ax1.plot(x, dat, color=(0, 0, 1))
dmean = np.nanmean(dat, axis=1)
ax1.plot(x, dmean, color=(1, 0, 0))
ax1.set_ylim([-1, 1])

ax2 = fig.add_subplot(2, 2, 2)
dim = 1  # y-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = corr_phase_xcom[:, dim, ft, ctcond, :]
ax2.plot(x, dat, color=(0, 0, 1))
dmean = np.nanmean(dat, axis=1)
ax2.plot(x, dmean, color=(1, 0, 0))
ax2.set_ylim([-1, 1])

ax3 = fig.add_subplot(2, 2, 3)
dim = 0  # x-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = gain_phase_xcom[:, dim, ft, ctcond, :]
ax3.plot(x, dat, color=(0, 0, 1))
dmean = np.nanmean(dat, axis=1)
ax3.plot(x, dmean, color=(1, 0, 0))
ax3.set_ylim([-2500, 1000])

ax4 = fig.add_subplot(2, 2, 4)
dim = 1  # y-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = gain_phase_xcom[:, dim, ft, ctcond, :]
ax4.plot(x, dat, color=(0, 0, 1))
dmean = np.nanmean(dat, axis=1)
ax4.plot(x, dmean, color=(1, 0, 0))
ax4.set_ylim([-2500, 1000])
plt.show()