%% Relation COM state and GRF (Part 2 tutorial)
%----------------------------------------------

clear all; close all; clc;
% settings
% datapath = 'C:\Users\mat950\Documents\Data\DataMoira_CSV\AnkleMoment';
datapath = 'C:\Users\mat950\Documents\Data\DataMoira_CSV\FP_Force_COM';
sVect = [3,4,7,8,14,18,21,22,24,25,31,36,41,44]; 
WalkCond = {'SteadyState_normal', 'SteadyState_slow'};
addpath(genpath(pwd));

% pre-allocate outcomes
nCond = length(WalkCond);
nSubj = length(sVect);
corr_phase_xcom = nan(50,2,2,nCond,nSubj);
gain_phase_xcom= nan(50,2,2,nCond,nSubj);
lag_xcom= nan(2,2,nCond,nSubj);
corr_phase_com= nan(50,2,2,nCond,nSubj);
gain_phase_com= nan(50,2,2,nCond,nSubj);
gain_phase_vcom= nan(50,2,2,nCond,nSubj);
lag_com= nan(2,2,nCond,nSubj);

cts = 0;
for s =sVect
    cts = cts + 1;
    ctcond = 0;
    for iCond = 1:length(WalkCond)
        ctcond = ctcond +1;
        % datapath
        filepath_data = fullfile(datapath, ['/S' num2str(s)],[WalkCond{iCond} '_data.csv']);
        filepath_event = fullfile(datapath, ['/S' num2str(s)],['/' WalkCond{iCond} '_event.csv']);
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
            [corr_phase_xcom(:,:,:,ctcond,cts),gain_phase_xcom(:,:,:,ctcond,cts),...
                lag_xcom(:,:,ctcond,cts),corr_phase_com(:,:,:,ctcond,cts), ...
                gain_phase_com(:,:,:,ctcond,cts),gain_phase_vcom(:,:,:,ctcond,cts),...
                lag_com(:,:,ctcond,cts)] = feedback_com_xcom(GRF,COM_filt,COMd,...
                xCOM,FootL_filt,FootR_filt,maxlag,Event.lhs,Event.rhs,time);
            
        end
    end
    disp(['Subject ' num2str(s) ' finished']);
end

%% plot results


figure();
subplot(2,2,1)
xVal = 51:100;
plot(xVal,squeeze(corr_phase_xcom(:,1,2,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
plot(xVal,nanmean(squeeze(corr_phase_xcom(:,1,2,1,:)),2),'Color',[0 0 0],'LineWidth',2); hold on;
set(gca,'YLim',[-1,1]);
subplot(2,2,2)
xVal = 51:100;
plot(xVal,squeeze(corr_phase_xcom(:,2,2,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
plot(xVal,nanmean(squeeze(corr_phase_xcom(:,2,2,1,:)),2),'Color',[0 0 0],'LineWidth',2); hold on;
set(gca,'YLim',[-1,1]);
subplot(2,2,3)
xVal = 51:100;
plot(xVal,squeeze(gain_phase_xcom(:,1,2,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
plot(xVal,nanmean(squeeze(gain_phase_xcom(:,1,2,1,:)),2),'Color',[0 0 0],'LineWidth',2); hold on;
set(gca,'YLim',[-2500,1000]);
subplot(2,2,4)
xVal = 51:100;
plot(xVal,squeeze(gain_phase_xcom(:,2,2,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',1.4); hold on;
plot(xVal,nanmean(squeeze(gain_phase_xcom(:,2,2,1,:)),2),'Color',[0 0 0],'LineWidth',2); hold on;
set(gca,'YLim',[-2500,1000]);
