function [Rsq, kp, kv, stats] = relate_com_anklemoment(t, com, foot, ankle_moment, ...
    events, varargin)
%relate_com_anklemoment Summary of this function goes here
%   Detailed explanation goes here

treadmill_velocity = NaN;
if ~isempty(varargin)
    treadmill_velocity = varargin{1};
end

% other inputs: feedback delay
fb_delay = 0.1; % feedback delay of 100 ms is default
x_stancephase = 0.05:0.05:0.95; % evaluate feedback at multiple instances in gait cycle
bool_nondim = false;
nhs_omit_end = 3; % omits first 3 gait cycles
boolplot = true;
BoolRemoveOutliers = true;

% get numerical derivatives
fs = 1./nanmean(diff(t));
comdot = calc_derivative(com,fs);
Foot_dot = calc_derivative(foot,fs);

% compute treadmill velocity if needed
if isnan(treadmill_velocity)
    % get all indices stance phase

    % average velocity foot during stance phase

end

% get COM position and velocity w.r.t foot
com_foot = com-foot;
comdot_foot = com+treadmill_velocity; % + because treadmill_velocity is absolute value

% nondim inputs if desired
if bool_nondim
    disp('todo')
end

% pre allocate matrices with dependent and independen variables
n_lhs = length(events.lhs);
COM_mat = nan(n_lhs, length(x_stancephase));
COMd_mat = nan(n_lhs, length(x_stancephase));
Tankle_mat = nan(n_lhs, length(x_stancephase));

% loop over all gait cyles to relate ankle moment and COM state
% loop over left heelstrikes
n_lhs = length(events.lhs);
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
    iSel =  X(:,1)>prctile(X(:,1),2) & X(:,1)<prctile(X(:,1),98) &...
        X(:,2)>prctile(X(:,2),2) & X(:,2)<prctile(X(:,2),98) &...
        Y>prctile(Y,2) & Y<prctile(Y,98);
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

% plot figure if desired
mk = 2;
if boolplot
    % standard figure
    h = figure('Color',[1 1 1],'Name','vis correlation');
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
    h = figure('Color',[1 1 1],'Name','vis correlation');
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
