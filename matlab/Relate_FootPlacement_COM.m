%% Relate foot placement and COM
%--------------------------------


clear all; close all; clc;
% path to experimental data
% info global coordinate system
%   x: medio-lateral
%   y: walking direction
%   z: vertical
datapath = '../ExampleData';

% settings
pred_samples = 1:50; % use ?
order = 2; % order of derivatives (1 = include velocity), (2= velocity and acc)
removeorigin = true; % 1= COM state w.r.t. contralateral foot
centerdata = true; % true = demean dependent and independent variables

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
    fs = 1./nanmean(diff(Dat.time)); % get sampling frequency
    COM = Dat.COMx;
    Rfoot = Dat.FootRx;
    Lfoot = Dat.FootLx;
    % convert events to index frame (instead of time)
    events.lhs = round(Event.lhs*fs + 1); % index starts at 1 in matlab so + 1
    events.rhs = round(Event.rhs*fs+ 1); % index starts at 1 in matlab so + 1
    events.rto = round(Event.rto*fs+ 1); % index starts at 1 in matlab so + 1
    events.lto = round(Event.lto*fs+ 1); % index starts at 1 in matlab so + 1
    % run foot placement code
    [OUT,intermediates]=foot_placement_model_function_step(COM,Rfoot, ...
        Lfoot,events,fs,pred_samples,order,removeorigin,centerdata,...
        'BoolPlot',true);
else
    if ~exist(filepath_data,'file')
        disp([filepath_data ' not on computer'])
    end
    if ~exist(filepath_event,'file')
        disp([filepath_event ' not on computer'])
    end
end
