%% Relate Ankle moment GRF
%-------------------------

clear all; close all; clc;
% settings
% path to experimental data
% info global coordinate system
%   x: medio-lateral
%   y: walking direction
%   z: vertical
datapath = 'C:\Users\mat950\Documents\Data\Moore2013_CSV\067';

% path to datafiles
filepath_data = fullfile(datapath,'UnpertPre_data.csv');
filepath_event = fullfile(datapath,'UnpertPre__event.csv');

% filepath_data = fullfile(datapath,'Pert_data.csv');
% filepath_event = fullfile(datapath,'Pert__event.csv');


if exist(filepath_data,'file') && exist(filepath_event,'file')
    % read data files
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event);
    t = Dat.time;

    % extract info (y is walking direction in this datast)
    COM = Dat.COM_x; % anterior-posterior COM motion
    FootL = Dat.LeftFoot_x;
    FootR = Dat.RightFoot_x;
    AnkleMoment = Dat.LeftAnklePlantarFlexionMoment; % plantarflexion/dorsiflexion

    % treadmill velocity in m/s
    treadmill_velocity = Dat.LeftBeltSpeed;% use nan if you want to compute this from marker coordinate
    
    % relate ankle moment to foot placement
    [stats] = relate_com_anklemoment(t, COM, FootL,AnkleMoment, Event, ...
        treadmill_velocity);
end
