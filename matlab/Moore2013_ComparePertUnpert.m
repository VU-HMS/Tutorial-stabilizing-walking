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

% path to datafiles
ConditionNames = {'UnpertPre','Pert','UnpertPost'};


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

% defualt plots
BoolDefaultPlots = false;

%% analysis
for ifile = 1:length(ConditionNames)
    % path to datafiles
    filepath_data = fullfile(datapath,[ConditionNames{ifile} '_data.csv']);
    filepath_event = fullfile(datapath,[ConditionNames{ifile} '__event.csv']);

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
        [output(ifile).corr_phase_xcom,output(ifile).gain_phase_xcom, ...
            output(ifile).lag_xcom,corr_phase_com, ...
            output(ifile).gain_phase_com, output(ifile).gain_phase_vcom,...
            output(ifile).lag_com] = feedback_com_xcom(GRF, ...
            COM_filt,FootL_filt,FootR_filt,events,fs,maxlag,L,...
            'BoolPlot', BoolDefaultPlots,'treadmill_velocity', treadmill_velocity_xz);

        % relate Foot placement to COM state
        [output(ifile).FootPlacement,output(ifile).FootPlacement_errors]= ...
            foot_placement_model_function_step(Dat.COM_z,Dat.RightFoot_z, ...
            Dat.LeftFoot_z,events,fs,pred_samples,order_derivative,removeorigin,centerdata,...
            'BoolPlot', BoolDefaultPlots,'treadmill_velocity', NaN); % nan veloctiy because medio-lat

        % relate ankle moment to COM state
        [output(ifile).Ta_Rsq, output(ifile).Ta_kp, output(ifile).Ta_kv,...
            output(ifile).resid, output(ifile).Ta_stats, output(ifile).Ta_settings] = ...
            relate_com_anklemoment(t, Dat.COM_x, Dat.LeftFoot_x,...
            Dat.LeftAnklePlantarFlexionMoment, Event, Anklemoment_delay, ...
            'treadmill_velocity',treadmill_velocity, 'BoolPlot',...
            BoolDefaultPlots, 'RemoveOutliers',true);
    else
        if ~exist(filepath_data,'file')
            disp([filepath_data ' not on computer'])
        end
        if ~exist(filepath_event,'file')
            disp([filepath_event ' not on computer'])
        end
    end
end


%% example plot compare three walking conditions

Colors_sel = [1 0 0; ...
    0 1 0;...
    0 0 1];

h_grf = figure('Name','GRF-COM analysis','Color',[1 1 1]);
h_foot = figure('Name','foot placement analysis','Color',[1 1 1]);
h_ankle = figure('Name','Ankle-COM analysis','Color',[1 1 1]);

for ifile = 1:length(ConditionNames)

    % select color
    Cs = Colors_sel(ifile,:);

    % plot relation COM / GRF
    %-----------------------------------
    figure(h_grf);
    subplot(2,2,1)
    xVal = 51:100;
    plot(xVal,squeeze(output(ifile).corr_phase_xcom(:,1,2)),'Color',Cs,'LineWidth',2); hold on;
    set(gca,'YLim',[-1,1]);
    ylabel('corr coef anterior-posterior')

    subplot(2,2,2)
    xVal = 51:100;
    plot(xVal,squeeze(output(ifile).corr_phase_xcom(:,2,2)),'Color',Cs,'LineWidth',2); hold on;
    set(gca,'YLim',[-1,1]);
    ylabel('corr coef medio-lateral')

    subplot(2,2,3)
    xVal = 51:100;
    plot(xVal,squeeze(output(ifile).gain_phase_xcom(:,1,2)),'Color',Cs,'LineWidth',2); hold on;
    xlabel('% gait cycle');
    ylabel('gain anterior-posterior')

    subplot(2,2,4)
    xVal = 51:100;
    plot(xVal,squeeze(output(ifile).gain_phase_xcom(:,2,2)),'Color',Cs,'LineWidth',2); hold on;
    xlabel('% gait cycle');
    ylabel('gain medio-lateral')
    
    if ifile == length(ConditionNames)
        for i =1:4
            subplot(2,2,i);
            set(gca,'box','off');
            set(gca,'LineWidth',1.6);
            set(gca,'FontSize',12);
        end
        legend(ConditionNames);
    end

    % plot foot placement analysis
    %-----------------------------------
    figure(h_foot);
    subplot(1,2,1);
    plot(output(ifile).FootPlacement.Combined_pct.data,'Color',Cs,'LineWidth',2);hold on;
    ylabel(output(ifile).FootPlacement.Combined_pct.titel);
    if ifile == length(ConditionNames)
        set(gca,'box','off');
        set(gca,'LineWidth',1.6);
        set(gca,'FontSize',12);
        xlabel('% gait cycle');
        legend(ConditionNames);
    end

    subplot(1,2,2);
    plot(nanmean(abs(output(ifile).FootPlacement_errors.error_combined),2),...
        'Color',Cs,'LineWidth',2);hold on;
    if ifile == length(ConditionNames)
        set(gca,'box','off');
        set(gca,'LineWidth',1.6);
        set(gca,'FontSize',12);
        xlabel('% gait cycle');
        ylabel('abs prediction error i.e. residual [m]')
        legend(ConditionNames);
    end

    % plot relation COM / ankle moment
    %-----------------------------------
    x_stancephase = output(1).Ta_settings.x_stancephase;
    figure(h_ankle);
    subplot(1,4,1)
    plot(x_stancephase, output(ifile).Ta_Rsq,'Color',Cs,'LineWidth',2); hold on;
    ylabel('Rsq')
    subplot(1,4,2)
    plot(x_stancephase, output(ifile).Ta_kp,'Color',Cs,'LineWidth',2); hold on;
    ylabel('position gain')
    subplot(1,4,3)
    plot(x_stancephase, output(ifile).Ta_kv,'Color',Cs,'LineWidth',2); hold on;
    ylabel('velocity gain')
    subplot(1,4,4)
    plot(x_stancephase, output(ifile).resid,'Color',Cs,'LineWidth',2); hold on;
    ylabel('abs residuals [Nm]');
    if ifile == length(ConditionNames)
        for i= 1:4
            subplot(1,4,i)
            set(gca,'box','off');
            set(gca,'LineWidth',1.6);
            set(gca,'FontSize',12);
            xlabel('% stance phase');
        end
        legend(ConditionNames);
    end
end
