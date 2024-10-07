%% Relate Ankle moment GRF
%-------------------------

%% path information
clear all; close all; clc;
% settings
% path to experimental data
% info global coordinate system
%   x: walking direction
%   y: vertical
%   z: medio-lateral
datapath = '../ExampleData/Moore2013';

% add path with functions
addpath(fullfile(pwd,'funcs'));

% path to datafiles -- unperturbed walking
% filepath_data = fullfile(datapath,'UnpertPre_data.csv');
% filepath_event = fullfile(datapath,'UnpertPre__event.csv');

% path to datafiles -- perturbed walking
filepath_data = fullfile(datapath,'Pert_data.csv');
filepath_event = fullfile(datapath,'Pert__event.csv');

%% Settings
% settings --- GRF =? COM
order = 2; % order lowpass filter
cutoff = 10; % cutoff lowpass filter
maxlag = 40; % maximal delay between COM state and GRF [number of frames]

% settings -- foot placement =? COM
pred_samples = 1:50; % use ?
order_derivative = 2; % order of derivatives (1 = include velocity), (2= velocity and acc)
removeorigin = true; % 1= COM state w.r.t. contralateral foot
centerdata = true; % true = demean dependent and independent variables

% settings --- AnkleMoment =? COM
Anklemoment_delay = 0.1; % 100 ms delay


%% analysis

if exist(filepath_data,'file') && exist(filepath_event,'file')
    % read data files
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event); % events
    t = Dat.time; % time vector
    L=nanmean(Dat.COM_y); % average height COM
    fs = 1./nanmean(diff(Dat.time)); % sampling frequency

    % convert events to index frame (instead of time): for foot placement
    % and GRF code
    events.lhs = round(Event.lhs*fs + 1); % index starts at 1 in matlab so + 1
    events.rhs = round(Event.rhs*fs+ 1); % index starts at 1 in matlab so + 1
    events.rto = round(Event.rto*fs+ 1); % index starts at 1 in matlab so + 1
    events.lto = round(Event.lto*fs+ 1); % index starts at 1 in matlab so + 1

    % get ground reaction forces and foot kinematics for GRF-COM code
    COM = [Dat.COM_x, Dat.COM_z];
    GRFL = [Dat.LeftGRF_x, Dat.LeftGRF_z];
    GRFR = [Dat.RightGRF_x, Dat.RightGRF_z];
    FootL = [Dat.LeftFoot_x Dat.LeftFoot_z]; % horizontal Foot position
    FootR = [Dat.RightFoot_x Dat.RightFoot_z];
    % get velocity of the treadmill
    treadmill_velocity = Dat.LeftBeltSpeed;% use nan if you want to compute this from marker coordinate

    % low pass filter the data
    COM_filt = LowpassFilterNan(COM,fs,order,cutoff);
    FootL_filt = LowpassFilterNan(FootL,fs,order,cutoff);
    FootR_filt = LowpassFilterNan(FootR,fs,order,cutoff);
    GRFR_filt = LowpassFilterNan(GRFR,fs,order,cutoff);
    GRFL_filt = LowpassFilterNan(GRFL,fs,order,cutoff);
    GRF = GRFR_filt + GRFL_filt; % combined GRF
    treadmill_velocity= LowpassFilterNan(treadmill_velocity, fs, order, cutoff);
    % add zero velocity in medio-lateral direciton
    treadmill_velocity_xz= [treadmill_velocity, zeros(length(t),1)];

    % relate GRF to COM state    
    [corr_phase_xcom,gain_phase_xcom, lag_xcom,corr_phase_com, ...
        gain_phase_com,gain_phase_vcom, lag_com] = feedback_com_xcom(GRF, ...
        COM_filt,FootL_filt,FootR_filt,events,fs,maxlag,L,...
        'BoolPlot', true,'treadmill_velocity', treadmill_velocity_xz);

    % relate Foot placement to COM state
    [OUT,intermediates]=foot_placement_model_function_step(Dat.COM_z,Dat.RightFoot_z, ...
        Dat.LeftFoot_z,events,fs,pred_samples,order_derivative,removeorigin,centerdata,...
        'BoolPlot', true,'treadmill_velocity', 0);

    % relate ankle moment to COM state    
    [OutTa.Rsq, OutTa.kp, OutTa.kv,OutTa.residuals, OutTa.stats] = ...
        relate_com_anklemoment(t, Dat.COM_x, Dat.LeftFoot_x,...
        Dat.LeftAnklePlantarFlexionMoment, Event, Anklemoment_delay, ...
        'treadmill_velocity',treadmill_velocity, 'BoolPlot', true, 'RemoveOutliers',true);
else
    if ~exist(filepath_data,'file')
        disp([filepath_data ' not on computer'])
    end
    if ~exist(filepath_event,'file')
        disp([filepath_event ' not on computer'])
    end
end
