%% Relate foot placement and COM
%--------------------------------


clear all; close all; clc;
% path to experimental data
datapath = 'C:\Users\mat950\Documents\Software\DataAnalysis\Tutorial-stabilizing-walking\ExampleData';

% settings
pred_samples = 1:50;
order = 2;
removeorigin = true;
centerdata = true;

% add subfolders
addpath(genpath(pwd));

% path to datafiles
filepath_data = fullfile(datapath,'S3','SteadyState_normal_data.csv');
filepath_event = fullfile(datapath,'S3','SteadyState_normal_event.csv');
if exist(filepath_data,'file') && exist(filepath_event,'file')
    % load data
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event);
    % restructure data
    fsopto = 1./nanmean(diff(Dat.time));
    COM = Dat.COMx;
    Rfoot = Dat.FootRx;
    Lfoot = Dat.FootLx;
    % convert events to index frame (instead of time)
    events.lhs = round(Event.lhs*fsopto + 1); % index starts at 1 in matlab so + 1
    events.rhs = round(Event.rhs*fsopto+ 1); % index starts at 1 in matlab so + 1
    events.rto = round(Event.rto*fsopto+ 1); % index starts at 1 in matlab so + 1
    events.lto = round(Event.lto*fsopto+ 1); % index starts at 1 in matlab so + 1
    % run foot placement code
    [OUT,intermediates]=foot_placement_model_function_step(COM,Rfoot,Lfoot,events,fsopto,pred_samples,order,removeorigin,centerdata);
end
% plot figure
figure();
plot(OUT.Combined_pct.data);
ylabel('variance explained');