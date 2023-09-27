import os

# use utilities
from utilities import *

# path information
datapath = "C:/Users/mat950/Documents/Software/DataAnalysis/Tutorial-stabilizing-walking/ExampleData"

# settings
centerdata = True
removeorigin = True
pred_samples = range(0, 50)
order = 2
boolplot = False


# loop over all subjects

# check if datafile exists
filepath_data = datapath + "/S3/SteadyState_normal_data.csv"
filepath_event = datapath + "/S3/SteadyState_normal_event.csv"
if os.path.exists(filepath_data) and os.path.exists(filepath_event):
    # read the data
    Dat = pd.read_csv(filepath_data)
    Event = pd.read_csv(filepath_event)

    # get inputs (mediolateral kinematics)
    t = np.array(Dat.time)
    Tankle = np.array(Dat.TAnkleLx)
    COM = np.array(Dat.COMx)
    FootL = np.array(Dat.FootLx)
    FootR = np.array(Dat.FootRx)

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
        boolplot)

# plot results
plt.figure()
plt.plot(pred_samples, output["Regression"]["Rsq"])
plt.xlabel('% single stance phase')
plt.ylabel("Variance explained")
#show figures
plt.show()