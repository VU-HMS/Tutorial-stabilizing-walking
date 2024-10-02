%% Relate Ankle moment GRF
%-------------------------

clear all; close all; clc;
% settings
% path to experimental data
% info global coordinate system
%   x: medio-lateral
%   y: walking direction
%   z: vertical
datapath = '../ExampleData';

% path to datafiles
filepath_data = fullfile(datapath,'S3','SteadyState_normal_data.csv');
filepath_event = fullfile(datapath,'S3','SteadyState_normal_event.csv');
if exist(filepath_data,'file') && exist(filepath_event,'file')
    % read data files
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event);
    t = Dat.time;
    L=nanmean(Dat.COMz); % average height COM

    % extract info (y is walking direction in this datast)
    COM = Dat.COMy; % anterior-posterior COM motion
    FootL = Dat.FootLy;
    FootR = Dat.FootRy;
    AnkleMoment = Dat.TAnkleLy; % plantarflexion/dorsiflexion

    % treadmill velocity in m/s
    treadmill_velocity = 1.1;% use nan if you want to compute this from marker coordinate
    
    % relate ankle moment to foot placement
    [stats] = relate_com_anklemoment(t, COM, FootL,AnkleMoment, Event, ...
        treadmill_velocity);
end