import os

# use utilities
from utilities import *


# had some issues with matplotlib in debug mode. this seems to solve all issues
import matplotlib
matplotlib.use('Qt5Agg') # interactive backend for matplotlib figures



# path information (point to the example data)
datapath = "../ExampleData"

# Create datatable for fit ?
BoolCreateTable = True

# init object to evaluate controller
# note that we work for no particular reason with an object oriented approach here.
# the idea behind it is that we store all datapoints in a large table and assign an
# id to each datapoint. This makes it easy to peform the regression on a specific subset
# of the data (e.g. for a specific phase in the gait cycle, of across subjects or for
# individual subjects).
myModel = reg_anklemoment_com(100000)

# time delay information
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
) # evaluate model at specific instances in the stance phase

# selected subjects
sVect = [3] # example data for subject 3 only

if BoolCreateTable:
    # check if datafile exists
    filepath_data = datapath + "/S3/SteadyState_normal_data.csv"
    filepath_event = datapath + "/S3/SteadyState_normal_event.csv"
    if os.path.exists(filepath_data) and os.path.exists(filepath_event):
        Dat = pd.read_csv(filepath_data) # read the datafile
        Event = pd.read_csv(filepath_event) # read the file with events
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
        Tanklefilt = ButterFilter_Low_NaNs(
            fs, Tankle, filteroder, filtercutoff
        )
        FootLfilt = ButterFilter_Low_NaNs(
            fs, FootL, filteroder, filtercutoff
        )
        # Numerical derivative for COM velocity
        COMdfilt = central_difference(t, COMfilt)

        # get treadmill speed
        FootdLfilt = central_difference(t, FootLfilt)
        trSpeed = np.nanmean(FootdLfilt[Dat.GRFLz>100]) # average velocity foot marker when foot in on the ground

        # nondim all inputs
        g = 9.81 # gravity
        m = np.mean(Dat.GRFLz + Dat.GRFRz) / g # mass subject
        L = np.mean(Dat.COMz) # height COM
        Tankle = Tanklefilt / (m * L * g)
        r_COM_Foot = (COMfilt - FootLfilt) / L # COM w.r.t. foot
        COMd = (COMdfilt-trSpeed) / np.sqrt(g * L)

        # loop over left heelstrikes
        n_lhs = len(Event.lhs)
        ctDelay = -1 # counter for phase in gait cycle
        for iDelay in ths_delay_vect:
            ctDelay = ctDelay + 1
            id = ctDelay # id connects the specific datapoints to a phase in the gait cycle
            ths_delay = iDelay
            for i in range(1, n_lhs - 3): # loop over all gait cycles
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
                # add datapoint to the regressionmodel object
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
    # load the data from a previously storded csv file
    myModel = reg_anklemoment_com(50000)
    myModel.importfromcsv("RegressionModel_DataMoira.csv")


# Run the regression and plot a figure
RsqAll_gaitPhase= np.zeros(shape=(len(ths_delay_vect)))
ctDelay = -1
BoolPlot = False # makes a plot for the regression model
for i in ths_delay_vect:
    ctDelay = ctDelay + 1
    id =ctDelay
    myModel.fitmodel(id, BoolPlot)  # fit model on all datapoints in this phase of the gait cycle
    RsqAll_gaitPhase[ctDelay ] = myModel.r2 # store the Rsq value
# plot the results
plt.figure()
plt.plot(ths_delay_vect*100,RsqAll_gaitPhase)
plt.xlabel('% stance phase')
plt.ylabel('Rsq')

# plot an example of the regression at a specific instance in the gait cycle
BoolPlot = True
id = 4 # at ths_delay_vect[id-1] in the gait cycle
myModel.fitmodel(id, BoolPlot)  # fit model on all datapoints in this phase of the gait cycle

# second test with function
fb_delay = 0.1
t = np.array(Dat.time)
Tankle = np.array(Dat.TAnkleLy)
COM = np.array(Dat.COMy)
FootL = np.array(Dat.FootLy)
relate_com_anklemoment(t, COM, FootL, Tankle, Event, fb_delay,
                       boolplot = True, treadmill_velocity= 1.1)

# show figures
plt.show()