import os

# use utilities
from utilities import *

# path information
datapath = "C:/Users/mat950/Documents/Software/DataAnalysis/Tutorial-stabilizing-walking/ExampleData"

# walking conditions
WalkCond = ["SteadyState_normal", "SteadyState_slow"]



# check if datafile exists
filepath_data = datapath + "/S3/SteadyState_normal_data.csv"
filepath_event = datapath + "/S3/SteadyState_normal_event.csv"
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
    [corr_phase_xcom, gain_phase_xcom, lag_xcom, corr_phase_com,
     gain_phase_com, gain_phase_vcom, lag_com] = feedback_com_xcom(
        GRF, COM_filt, COMd_filt, xCOM, FootL_filt, FootR_filt, maxlag, t_lhs, t_rhs, time
    )

fig = plt.figure()

x = np.linspace(51, 100, 50)
ax1 = fig.add_subplot(2, 2, 1)
dim = 0  # x-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = corr_phase_xcom[:, dim, ft]
ax1.plot(x, dat, color=(0, 0, 1))
ax1.set_ylim([-1, 1])
#ax1.set_xlabel('% cycle')
ax1.set_ylabel('R2')

ax2 = fig.add_subplot(2, 2, 2)
dim = 1  # y-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = corr_phase_xcom[:, dim, ft]
ax2.plot(x, dat, color=(0, 0, 1))
ax2.set_ylim([-1, 1])
#ax2.set_xlabel('% cycle')
#ax2.set_ylabel('R2')

ax3 = fig.add_subplot(2, 2, 3)
dim = 0  # x-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = gain_phase_xcom[:, dim, ft]
ax3.plot(x, dat, color=(0, 0, 1))
ax3.set_ylim([-2500, 1000])
ax3.set_xlabel('% cycle')
ax3.set_ylabel('gain')

ax4 = fig.add_subplot(2, 2, 4)
dim = 1  # y-direction
ft = 1  # right foot
ctcond = 0  # normal walking
dat = gain_phase_xcom[:, dim, ft]
ax4.plot(x, dat, color=(0, 0, 1))
ax4.set_ylim([-2500, 1000])
ax4.set_xlabel('% cycle')
#ax4.set_ylabel('gain')
plt.show()