# Tutorial-stabilizing-walking

This repository contains the matlab and python implementation of the methods presented in "Assessment of stabilizing feedback control of walking, a tutorial" (https://doi.org/10.1101/2023.09.28.559516). In this tutorial paper we analyse the relation between COM state and (1) ground reaction force, (2) foot placement and (3) ankle moments during walking.

*Note for users:* the python code has only been tested on the example data and is probably less robust

## Control of horizontal ground reaction forces

You can analyse the between relation center of mass kinematics and horizontal ground reaction forces and using the function **feedback_com_xcom** (matlab: funcs/feedback_com_xcom.m, python: utilities.py/feedback_com_xcom). The scripts Relate_GRF_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.

matlab example:

```matlab
[corr_xcom,gain_xcom, lag_xcom,corr_com, gaincom,gain_vcom, lag_com] = ...
feedback_com_xcom(GRF, COM, FootL, FootR, events, fs, maxlag, ...
L_pendulum, 'BoolPlot', true, 'treadmill_velocity',treadmill_velocity); 
```

open the function **feedback_com_xcom.m** for more details on the  required and optional input arguments and outputs.

## Control of foot placement

You can analyse the relation between center of mass kinematics and foot placement using the function **foot_placement_model_function_step** (matlab: funcs/foot_placement_model_function_step.m, python: utilities.py/FootPlacement_COMstate). The scripts Relate_FootPlacement_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.

matlab example:

```matlab
[OUT,intermediates]=foot_placement_model_function_step(COM,Rfoot,...
Lfoot,events,fs,pred_samples,order,removeorigin,centerdata,...
'BoolPlot',true);
```

open the function **foot_placement_model_function_step** for more details on the  required and optional input arguments and outputs.

## Control of ankle moment

You can analyse the relation between center of mass kinematics and ankle moment using the function relate_com_anklemoment (matlab: funcs/relate_com_anklemoment.m, python: class reg_anklemoment_com). The scripts Relate_AnkleMoment_COM.m and AnkleMoment_Unpert.py contain examples on how to run this analysis in matlab or python.

matlab example:

```matlab
[Rsq, Kp, Kv,  residuals, stats] = relate_com_anklemoment(t,...
 COM, FootL,AnkleMoment, Event, delay, 'treadmill_velocity',...
 treadmill_velocity, 'BoolPlot',true);
```

## Datasets

*unperturbed walking:* The examples use data of one subject (ExampleData/S3) of the study: [Dataset and analyses for: Active foot placement control ensures stable gait: Effect of constraints on foot placement and ankle moments](https://doi.org/10.5281/zenodo.4229851). This analysis can be applied to all subjects after downloading the data.

*perturbed and unperturbed walking*: the examples GeneralExample_Moore2013.m and Moore2013_ComparePertUnpert.m use the data from: Moore JK, Hnat SK, van den Bogert AJ. 2015. An elaborate data set on human gait and the effect of mechanical perturbations. PeerJ 3:e918 [An elaborate data set on human gait and the effect of mechanical perturbations [PeerJ]](https://doi.org/10.7717/peerj.918). Example data of one subject walking (ExampleData/Moore2013) is used in these example. The analysis can be applied to all subjects after dowloading the data from the online repository.

## Python examples

python example: Control of horizontal ground reaction forces

```python
maxlag = 40 # maximal lag is 40% of gait cycle
[corr_phase_xcom, gain_phase_xcom, lag_xcom, corr_phase_com,
 gain_phase_com, gain_phase_vcom, lag_com] = feedback_com_xcom(
    GRF, COM, COMd, xCOM, FootL,
    FootR, maxlag, t_lhs, t_rhs, time)
```

python example: Control of foot placement

```python
# settings
centerdata = True
removeorigin = True
pred_samples = range(0, 50)
order = 2

# run the foot placement model on this dataset
output = FootPlacement_COMstate(t,COM,FootL,FootR,Event,
                                centerdata,removeorigin,
                                pred_samples,order,boolplot)
```

python example control of ankle moment: You can analyse the relation between center of mass kinematics and ankle using the class reg_anklemoment_com (utilities.py/reg_anklemoment_com). The script AnkleMoment_Unpert.py contains an example on how to run this analysis in python

```python

# create object for regression model

myModel = reg_anklemoment_com(100000)

# example on how to add datapoint to object

id = 1 # identified for datapoint that enables running regression on subset of data
myModel.adddatapoint(COM, COMd, Tankle, id)

# fit regression

idsel = [1]
BoolPlot = True # make figure with datapoints and regression lines
myModel.fitmodel(idsel, BoolPlot)  # fit model on all with ids in idsel
VarExpl = myModel.r2 # get variance explained by regression model
