function [COM_filt] = LowpassFilterNan(COM,fs,order,cutoff)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% detect nans
iNan = isnan(COM(:,1));
iNoNaN = ~isnan(COM(:,1));

% interpolate the temporary remove nans
IndexAll = 1:length(COM);
IndexNoNaN = IndexAll(iNoNaN);
COM_int = interp1(IndexNoNaN, COM(iNoNaN,:),IndexAll, 'spline');

% filter data without nans
[a,b] = butter(order,cutoff./(0.5*fs),'low');
COM_filt = filtfilt(a,b,COM_int);
COM_filt(iNan,:) = NaN;

end