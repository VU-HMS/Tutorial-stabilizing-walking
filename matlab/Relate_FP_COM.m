%% Relate foot placement and COM
%--------------------------------


clear all; close all; clc;
% settings
% datapath = 'C:\Users\mat950\Documents\Data\DataMoira_CSV\AnkleMoment';
datapath = 'C:\Users\mat950\Documents\Data\DataMoira_CSV\FP_Force_COM';
% sVect = [3,4,7,8,14,18,21,22,24,25,31,36,41,44];
sVect = 3;
WalkCond = {'SteadyState_normal', 'SteadyState_slow'};
addpath(genpath(pwd));

% pre-allocate outcomes
Rsq_RightLeg = nan(50,length(sVect),2);
Rsq_LeftLeg = nan(50,length(sVect),2);
Rsq_Legs = nan(50,length(sVect),2);

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
            % load data
            Dat = readtable(filepath_data);
            Event = readtable(filepath_event);
            % restructure data
            fsopto = 1./nanmean(diff(Dat.time));
            pred_samples = 1:50;
            order = 2;
            removeorigin = true;
            centerdata = true;
            COM = Dat.COMx;
            Rfoot = Dat.FootRx;
            Lfoot = Dat.FootLx;
            events.lhs = round(Event.lhs*50 + 1); % index starts at 1 in matlab so + 1
            events.rhs = round(Event.rhs*50+ 1); % index starts at 1 in matlab so + 1
            events.rto = round(Event.rto*50+ 1); % index starts at 1 in matlab so + 1
            events.lto = round(Event.lto*50+ 1); % index starts at 1 in matlab so + 1
            % run foot placement code
            [OUT,intermediates]=foot_placement_model_function_step(COM,Rfoot,Lfoot,events,fsopto,pred_samples,order,removeorigin,centerdata);
            Rsq_RightLeg(:,cts,ctcond) = OUT.Right_pct.data;
            Rsq_LeftLeg(:,cts,ctcond) = OUT.Left_pct.data;
            Rsq_Legs(:,cts,ctcond) = OUT.Combined_pct.data;
        end
    end
    disp(['finished with subject ' num2str(s)]);
end

%%
figure();
subplot(1,2,1)
plot(Rsq_RightLeg(:,:,1));
subplot(1,2,2)
plot(Rsq_RightLeg(:,:,2));

figure();
subplot(1,2,1)
plot(Rsq_LeftLeg(:,:,1));
subplot(1,2,2)
plot(Rsq_LeftLeg(:,:,2));

figure();
subplot(1,2,1)
plot(Rsq_Legs(:,:,1));
subplot(1,2,2)
plot(Rsq_Legs(:,:,2));