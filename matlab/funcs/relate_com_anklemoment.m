function [Rsq, kp, kv, stats] = relate_com_anklemoment(t, com, foot, ankle_moment, ...
    events,fb_delay, varargin)
%relate_com_anklemoment computes correlation between delayed com position
%and velocity (w.r.t. stance foot) and the ankle joint moment (feel free to
% use another signal).
%   input arguments:
%       (1): t: time vector
%       (2): com: array with com position
%       (3): foot: array with position foot
%       (4): ankle_moment: array with ankle moment
%       (5): event: structure with events in time [s]
%               lhs = left heel strike
%               rhs = right heel strike
%               lto = left toe off
%               rto = right toe off
%       (6) fb_delay: feedback delay
%       (7) optional arguments:
%       - treadmill_velocity: velocity of treadmill (can be a scalar of
%                             array). If treadmill velocity is specified it
%                             is used to compute the velocity of the COM 
%                              w.r.t. foot (instead off numericaldervivative
%                              foot velocity).
%       - BoolPlot: creates default plot out outcomes
%       - stancephase: array ranging between 0 and 1 with descrete time for
%                      the feedback analysis.
%       - nhs_omit: omits the last n gait cycles (e.g. end of experiment)
%       - RemoveOutliers: if true datapoints outside 2 and 98th percentile 
%                         are removed fromt he analysis
%       


%% handle optional  inputs
% default values of optional input arguments
x_stancephase = 0.05:0.05:0.95; % evaluate feedback at multiple instances in gait cycle
bool_nondim = false;
nhs_omit_end = 0; % omits first 3 gait cycles
boolplot = true;
RemoveOutliers = true;
treadmill_velocity = NaN;

% Parse optional name-value pairs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'treadmill_velocity'
            treadmill_velocity = varargin{i+1};
        case 'BoolPlot'
            BoolPlot = varargin{i+1};
        case 'stancephase'
            x_stancephase = varargin{i+1};
        case 'nhs_omit'
            nhs_omit = varargin{i+1};
        case 'RemoveOutliers'
            RemoveOutliers = varargin{i+1};
        otherwise
            error('Unknown parameter %s', varargin{i});
    end
end

%% compute position and velocity com w.r.t. stance foot
% get numerical derivatives
fs = 1./nanmean(diff(t));
comdot = calc_derivative(com,fs);
Foot_dot = calc_derivative(foot,fs);

% number of left heelstrikes
n_lhs = length(events.lhs);

% compute treadmill velocity if needed
if isnan(treadmill_velocity)
    % get all indices single stance phase
    istance= zeros(length(t),1);
    for i=1:length(events.rto)
        if ~isnan(events.rto(i))
            % get first rhs after rto 
            iabove = events.rhs > events.rto(i);
            tend = events.rhs(find(iabove, 1));
            % indices single stance phase left leg
            isel = t>events.rto(i) & t<tend ;
            istance(isel) = 1;
        end
    end
    % computes average velocity treadmill during single stance phase
    % be aware that this assumes a constant treadmill velocity
    treadmill_velocity = abs(nanmean(Foot_dot(istance == 1)));
    % abs because we use absolute value of velocity (otherwise the input
    % arguments should typically be a negative number).
end

% get COM position and velocity w.r.t foot
com_foot = com-foot;
comdot_foot = com+treadmill_velocity; % + because treadmill_velocity is absolute value

% nondim inputs if desired
if bool_nondim
    disp('todo')
end

%% create table with delayed COM state and ankle moment at descrete time points
% pre allocate matrices with dependent and independen variables
n_lhs = length(events.lhs);
COM_mat = nan(n_lhs, length(x_stancephase));
COMd_mat = nan(n_lhs, length(x_stancephase));
Tankle_mat = nan(n_lhs, length(x_stancephase));

% loop over all gait cyles to relate ankle moment and COM state
% loop over left heelstrikes
ctDelay = 0; % counter for phase in gait cycle
for iDelay = x_stancephase
    ctDelay = ctDelay + 1;
    ths_delay = iDelay;
    for i = 1:(n_lhs - nhs_omit_end) % loop over all gait cycles
        % find toe-off
        iabove = events.lto > events.lhs(i);
        dt_stance = events.lto(find(iabove, 1)) - events.lhs(i);
        if ~isnan(events.lhs(i))
            % get COM state information
            iabove = t > (events.lhs(i) + ths_delay * dt_stance);
            ix = find(iabove, 1);
            COM_mat(i, ctDelay) = com_foot(ix);
            COMd_mat(i, ctDelay) = comdot_foot(ix);

            % get Torque information
            iabove = t > (events.lhs(i) + ths_delay * dt_stance + fb_delay);
            iy = find(iabove, 1);
            Tankle_mat(i, ctDelay) = ankle_moment(iy);
        end
    end
end

%% correlation
% compute regression
for i=1:length(x_stancephase)
    % get dependent and independent variables
    X = [COM_mat(:,i), COMd_mat(:,i)];
    Y = Tankle_mat(:,i);
    % remove NaNs
    inan = isnan(X(:,1)) | isnan(X(:,2)) | isnan(Y);
    X(inan,:) = [];
    Y(inan,:) = [];
    % ony include data in 5-95 precentile
    if RemoveOutliers
        iSel =  X(:,1)>prctile(X(:,1),2) & X(:,1)<prctile(X(:,1),98) &...
            X(:,2)>prctile(X(:,2),2) & X(:,2)<prctile(X(:,2),98) &...
            Y>prctile(Y,2) & Y<prctile(Y,98);
    end
    X = X(iSel,:);
    Y = Y(iSel);
    % run statistics
    stats(i) = regstats(Y,X,'linear',{'beta','covb','yhat','r',...
        'mse','rsquare','adjrsquare','tstat','fstat'});
    % store inputs
    inputs(i).X = X;
    inputs(i).Y= Y;

end

% create output structure
Rsq = [stats().rsquare];
gains = [stats().beta];
kp = gains(2,:);
kv = gains(3,:);

%% default figure

% plot figure if desired
mk = 2;
if boolplot
    % standard figure
    h = figure('Color',[1 1 1],'Name','Relate Moment-COM: outputs');
    subplot(1,3,1)
    plot(x_stancephase, Rsq,'Color',[0 0 0],'LineWidth',2);
    ylabel('Rsq')
    subplot(1,3,2)
    plot(x_stancephase, kp,'Color',[0 0 0],'LineWidth',2);
    ylabel('position gain')
    subplot(1,3,3)
    plot(x_stancephase, kv,'Color',[0 0 0],'LineWidth',2);
    ylabel('velocity gain')
    for i= 1:3
        subplot(1,3,i)
        set(gca,'box','off');
        set(gca,'LineWidth',1.6);
        set(gca,'FontSize',12);
        xlabel('% stance phase');
    end
    % a bit more adventurous figure to explore ankle moment and variance
    % explained
    y_scale = (max(ankle_moment) - min(ankle_moment))*7;
    h = figure('Color',[1 1 1],'Name','Relate Moment-COM: interpret relation');
    for i=1:length(x_stancephase)
        % get dependent and independent variables
        Y = inputs(i).Y;
        x = x_stancephase(i);
        plot(x + stats(i).yhat./y_scale, Y./y_scale,'ok','Color',[0 0 0],...
            'MarkerFaceColor',[0 0 0],'MarkerSize',mk); hold on;
        plot(x + [min(stats(i).yhat),max(stats(i).yhat)]./y_scale,...
            [min(Y./y_scale) max(Y./y_scale)],'Color',[0.6 0.6 0.6],'LineWidth',1);
        xlabel('% stance phase (and var recon norm ankle moment [])');
        ylabel('norm ankle moment []')
    end
    set(gca,'box','off');
    set(gca,'LineWidth',1.6);
    set(gca,'FontSize',12);
end

end
