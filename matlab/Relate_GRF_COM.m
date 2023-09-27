%% Relation COM state and GRF (Part 2 tutorial)
%----------------------------------------------

clear all; close all; clc;
% settings
% path to experimental data
datapath = 'C:\Users\mat950\Documents\Software\DataAnalysis\Tutorial-stabilizing-walking\ExampleData';

% path to datafiles
filepath_data = fullfile(datapath,'S3','SteadyState_normal_data.csv');
filepath_event = fullfile(datapath,'S3','SteadyState_normal_event.csv');
if exist(filepath_data,'file') && exist(filepath_event,'file')
    Dat = readtable(filepath_data);
    Event = readtable(filepath_event);

    % extract info
    COM = [Dat.COMx Dat.COMy];
    time = Dat.time;
    FootL = [Dat.FootLx Dat.FootLy];
    FootR = [Dat.FootRx Dat.FootRy];
    GRFR = [Dat.GRFRx Dat.GRFRy];
    GRFL = [Dat.GRFLx Dat.GRFLy];

    % low pass filter the data
    fs = 1./nanmean(diff(time));
    order = 2;
    cutoff = 10;
    COM_filt = LowpassFilterNan(COM,fs,order,cutoff);
    FootL_filt = LowpassFilterNan(FootL,fs,order,cutoff);
    FootR_filt = LowpassFilterNan(FootR,fs,order,cutoff);
    GRFR_filt = LowpassFilterNan(GRFR,fs,order,cutoff);
    GRFL_filt = LowpassFilterNan(GRFL,fs,order,cutoff);
    GRF = GRFR_filt + GRFL_filt;

    % compute COM velocity and xCOM position
    COMd = NumericalDerivative_CentralDiff(time, COM_filt);
    L = nanmean(Dat.COMz);
    xCOM = COM_filt + COMd./ sqrt(9.81/L);
    maxlag = 40;
    [corr_phase_xcom,gain_phase_xcom, lag_xcom,corr_phase_com, ...
        gain_phase_com,gain_phase_vcom, lag_com] = feedback_com_xcom(GRF,COM_filt,COMd,...
        xCOM,FootL_filt,FootR_filt,maxlag,Event.lhs,Event.rhs,time);
end


%% plot results


figure();
subplot(2,2,1)
xVal = 51:100;
plot(xVal,squeeze(corr_phase_xcom(:,1,2)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
set(gca,'YLim',[-1,1]);
subplot(2,2,2)
xVal = 51:100;
plot(xVal,squeeze(corr_phase_xcom(:,2,2)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
set(gca,'YLim',[-1,1]);
subplot(2,2,3)
xVal = 51:100;
plot(xVal,squeeze(gain_phase_xcom(:,1,2)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
set(gca,'YLim',[-2500,1000]);
subplot(2,2,4)
xVal = 51:100;
plot(xVal,squeeze(gain_phase_xcom(:,2,2)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
set(gca,'YLim',[-2500,1000]);
