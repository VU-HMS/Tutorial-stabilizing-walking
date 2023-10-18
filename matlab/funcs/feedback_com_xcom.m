function [corr_phase_xcom,gain_phase_xcom,lag_xcom,corr_phase_com,gain_phase_com,...
    gain_phase_vcom,lag_com] = feedback_com_xcom(force,com,vcom,xcom,lfoot,rfoot,...
    maxlag,l_HeelStrike,r_HeelStrike,time)
% feedback_com_xcom computed the relation between the ground reaction
% forces and the extrapolated center of mass information as a function of
% the gait cycle for a range of delays  between both signals
% Inputs:
%   (1) force: nx3 vector with ground reaction forces
%   (2) com: nx3 vector with center of mass position
%   (3) vcom: nx3 vector with center of mass velocity
%   (4) xcom: nx3 vector with extrapolated center of mass position
%   (5) lfoot: nx3 vector with position of the left foot
%   (6) lfoot: nx3 vector with position of the right foot
%   (7) max: maximal lag in the system (expressed as % gait cycle, so not in s)
%   (8) l_Heelstrike: vector with timing of the left heelstrikes
%   (9) r_Heelstrike: vector with timing of the right heelstrikes
%   (10) time: time vector


n_strides = length(l_HeelStrike) - 1;
stride_f = nanmean(diff(l_HeelStrike));

% pre-allocate outcomes
corr_phase_xcom = nan(50,2,2);
gain_phase_xcom = nan(50,2,2);
corr_phase_com = nan(50,2,2);
gain_phase_com = nan(50,2,2);
gain_phase_vcom = nan(50,2,2);
lag_xcom = nan(2,2);
lag_com = nan(2,2);


for dim = 1:2 % dim = 1 for ML; dim = 2 for AP
    dimdata  = [force(:,dim),com(:,dim),vcom(:,dim),xcom(:,dim),lfoot(:,dim),rfoot(:,dim)];
    for ft = 1:2
        % calculate time_normalized forces, xcom, etc
        nLeftHS = length(l_HeelStrike);
        nRightHS = length(r_HeelStrike);
        maxHs = max([nLeftHS nRightHS]);
        ndata2 = nan(100,maxHs-1,6);
        for var = 1:6
            if ft == 1
                %tmp = normalizetimebase(dimdata(:,var),l_HeelStrike);
                tmp = NormCycle(time,l_HeelStrike,dimdata(:,var));
            else
                tmp = NormCycle(time,r_HeelStrike,dimdata(:,var));
            end
            ncycl = length(tmp(1,:));
            ndata2(:,1:ncycl,var)  = tmp;
        end
        % reference com and xcom to foot at start of cycle
        ndata2(:,:,2)        =  ndata2(:,:,2) - nanmean(ndata2(10:40,:,4+ft));
        ndata2(:,:,4)        =  ndata2(:,:,4) - nanmean(ndata2(10:40,:,4+ft));
        % association between forces in stride with xcom at phase lead j
        corr_force_xcom = nan(50,maxlag);
        gain_force_xcom = nan(50,maxlag);
        gain_force_com = nan(50,maxlag);
        gain_force_vcom = nan(50,maxlag);
        corr_force_com = nan(50,maxlag);
        for j = 1:maxlag  % vary delay from 1 to maxlag% of gait cycle
            for i = 51:100 %phase
                mf = ndata2(i,:,1);
                mx = ndata2(i-j,:,4);
                mc = ndata2(i-j,:,2);
                mv = ndata2(i-j,:,3);
                iNan = isnan(mf) | isnan(mx);
                % if j == maxlag && i == 100
                %     Redisp('test');
                % end
                mf(iNan) = [];
                mx(iNan) = [];
                mc(iNan) = [];
                mv(iNan) = [];

                corr_force_xcom(i-50,j) = corr(mf',mx');
                p_temp = polyfit(mx,mf,1); 
                gain_force_xcom(i-50,j) = p_temp(1);

                stats                    = regstats(mf',[mc',mv'],'linear',{'beta','rsquare'});
                gain_force_com(i-50,j)   = stats.beta(2);
                gain_force_vcom(i-50,j)  = stats.beta(3);
                corr_force_com(i-50,j)   = sqrt(stats.rsquare);

            end
        end
        % find lag with strongest neg. correlation overall
        test            = corr_force_xcom;
        test(test>0)    = NaN;
        ztest           = tanh(nanmean(atanh(test),1));
        [~,lag_xcom(dim,ft)] = min(ztest);
        corr_phase_xcom(:,dim,ft)   = corr_force_xcom(:,lag_xcom(dim,ft));
        gain_phase_xcom(:,dim,ft)   = gain_force_xcom(:,lag_xcom(dim,ft));

        % find lag with most negative weighted gain
        beta_vel = squeeze(gain_force_vcom(:,:));
        beta_pos = squeeze(gain_force_com(:,:));
        beta_weighted               = beta_vel.*(2*pi*stride_f) + beta_pos;
        [~,lag_com(dim,ft)]         = min(min(beta_weighted));  % this needs to be improved it seems
%         [~,lag_com(dim,ft)]         = max(max(corr_force_com));  % this needs to be improved it seems
        corr_phase_com(:,dim,ft)    = corr_force_com(:,lag_com(dim,ft));
        gain_phase_com(:,dim,ft)    = gain_force_com(:,lag_com(dim,ft));
        gain_phase_vcom(:,dim,ft)   = gain_force_vcom(:,lag_com(dim,ft));

    end
end
