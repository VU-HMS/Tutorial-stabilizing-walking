# Tutorial-stabilizing-walking

This repository contains the matlab and python implementation of the methods presented in "Assessment of stabilizing feedback control of walking, a tutorial" [add doi]. In this tutorial paper we analyse the relation between COM state and (1) ground reaction force, (2) foot placement and (3) ankle moments during walking.

## Control of horizontal ground reaction forces

You can analyse the between relation center of mass kinematics and horizontal ground reaction forces and using the function **feedback_com_xcom** (matlab: funcs/feedback_com_xcom.m, python: utilities.py/feedback_com_xcom). The scripts Relate_FP_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.

## Control of foot placement

You can analyse the relation between center of mass kinematics and foot placement using the function **foot_placement_model_function_step** (matlab: funcs/foot_placement_model_function_step.m, python: utilities.py/FootPlacement_COMstate). The scripts Relate_FP_COM.m and FootPlacement_Unpert.py contain examples on how to run this analysis in matlab or python.

## Control of ankle moment

You can analyse the relation between center of mass kinematics and ankle using the class reg_anklemoment_com (utilities.py/reg_anklemoment_com) . The script AnkleMoment_Unpert.py contains an example on how to run this analysis in python



## Datasets

The examples work with the data of one subject (ExampleData) but can be applied to all subjects after downloading the data from: https://doi.org/10.5281/zenodo.4229851 (Dataset and analyses for: Active foot placement control ensures stable gait: Effect of constraints on foot placement and ankle moments)




