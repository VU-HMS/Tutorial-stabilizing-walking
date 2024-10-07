%% Relation COM state and GRF (Part 2 tutorial)
%----------------------------------------------

clear all; close all; clc;
% settings
% path to experimental data
% info global coordinate system
%   x: medio-lateral
%   y: walking direction
%   z: vertical
datapath = '../ExampleData';

% add path with functions
addpath(fullfile(pwd,'funcs'));

% path to datafiles
filepath_data = fullfile(datapath,'S3','SteadyState_normal_data.csv');
filepath_event = fullfile(datapath,'S3','SteadyState_normal_event.csv');
if exist(filepath_data,'file') && exist(filepath_event,'file')
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event);
    fs = 1./nanmean(diff(Dat.time)); % sampling frequency
    L=nanmean(Dat.COMz); % average height COM

    % extract info (y is walking direction in this datast)
    COM = [Dat.COMx Dat.COMy]; % horizontal COM motion
    FootL = [Dat.FootLx Dat.FootLy]; % horizontal Foot position
    FootR = [Dat.FootRx Dat.FootRy];
    GRFR = [Dat.GRFRx Dat.GRFRy];
    GRFL = [Dat.GRFLx Dat.GRFLy];
    % note that GRF was already resampled to motion capture framerate

    % convert events to index frame (instead of time)
    events.lhs = round(Event.lhs*fs + 1); % index starts at 1 in matlab so + 1
    events.rhs = round(Event.rhs*fs+ 1); % index starts at 1 in matlab so + 1
    events.rto = round(Event.rto*fs+ 1); % index starts at 1 in matlab so + 1
    events.lto = round(Event.lto*fs+ 1); % index starts at 1 in matlab so + 1

    % low pass filter the data
    order = 2;
    cutoff = 10;
    COM_filt = LowpassFilterNan(COM,fs,order,cutoff);
    FootL_filt = LowpassFilterNan(FootL,fs,order,cutoff);
    FootR_filt = LowpassFilterNan(FootR,fs,order,cutoff);
    GRFR_filt = LowpassFilterNan(GRFR,fs,order,cutoff);
    GRFL_filt = LowpassFilterNan(GRFL,fs,order,cutoff);
    GRF = GRFR_filt + GRFL_filt; % combined GRF

    % treadmill_velocity
    treadmill_velocity = [0 1.1]; % 0 in medio-lateral, 1.1 in anterior_posterior

    % relate GRF to COM state
    maxlag = 40; % maximal delay between COM state and GRF [number of frames]
    [corr_phase_xcom,gain_phase_xcom, lag_xcom,corr_phase_com, ...
        gain_phase_com,gain_phase_vcom, lag_com] = feedback_com_xcom(GRF, ...
        COM_filt,FootL_filt,FootR_filt,events,fs,maxlag,L, 'BoolPlot', true,...
        'treadmill_velocity',treadmill_velocity);
else
    if ~exist(filepath_data,'file')
        disp([filepath_data ' not on computer'])
    end
    if ~exist(filepath_event,'file')
        disp([filepath_event ' not on computer'])
    end
end
