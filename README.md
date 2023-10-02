# Tutorial-stabilizing-walking

This repository contains the matlab and python implementation of the methods presented in "Assessment of stabilizing feedback control of walking, a tutorial" (https://doi.org/10.1101/2023.09.28.559516). In this tutorial paper we analyse the relation between COM state and (1) ground reaction force, (2) foot placement and (3) ankle moments during walking.

## Control of horizontal ground reaction forces

You can analyse the between relation center of mass kinematics and horizontal ground reaction forces and using the function **feedback_com_xcom** (matlab: funcs/feedback_com_xcom.m, python: utilities.py/feedback_com_xcom). The scripts Relate_FP_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.



```python
maxlag = 40 # maximal lag is 40% of gait cycle
[corr_phase_xcom, gain_phase_xcom, lag_xcom, corr_phase_com,
 gain_phase_com, gain_phase_vcom, lag_com] = feedback_com_xcom(
    GRF, COM, COMd, xCOM, FootL,
    FootR, maxlag, t_lhs, t_rhs, time)
```



## Control of foot placement

You can analyse the relation between center of mass kinematics and foot placement using the function **foot_placement_model_function_step** (matlab: funcs/foot_placement_model_function_step.m, python: utilities.py/FootPlacement_COMstate). The scripts Relate_FP_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.



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



## Control of ankle moment

You can analyse the relation between center of mass kinematics and ankle using the class reg_anklemoment_com (utilities.py/reg_anklemoment_com) . The script AnkleMoment_Unpert.py contains an example on how to run this analysis in python



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

```



## Datasets

The examples work with the data of one subject (ExampleData) but can be applied to all subjects after downloading the data from: https://doi.org/10.5281/zenodo.4229851 (Dataset and analyses for: Active foot placement control ensures stable gait: Effect of constraints on foot placement and ankle moments)




