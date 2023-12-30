function  [Cycle,TimeGain]=NormalizeTimeBase(signal,trigind)

% NormalizeTimeBase: calcultes a 0-100% timebase, ensemble-averages cyclic signals
% [Cycle,TimeGain]=NormalizeTimeBase(signal,trigind)
%
% Input : Signal: any one-dimensional array
%         trigind : array of indices, default: [1 length(Signal)]
%                   should increase monotoneusly
% Process: calculates new points based on a 0-100% time base
%          by spline interpolation for each time interval
% Output:  if length(trigind)=2: Cycle [100 1]
%          if length(trigind)>2: Cycle [100 Ncycles+2]
%             Ncyles=length(trigind)-1,
%             Ncycles+1: mean signal per point, i.e. ensemble averaged
%             Ncycles+2: stand.dev ensemble averaged points
%          TimeGain: (average) amplification/reduction of time-axis (i.e. 100/(samples/cycle))
%
% WARNING user should be aware of information loss in case of excessive downsampling

% AUTHOR(S) AND VERSION-HISTORY
% Ver 1.2 April 2003 (Jaap Harlaar VUmc Amsterdam) adapted from some version


if nargin < 2, trigind=[1 length(signal)]; end
if nargin < 1, return, end

% FFE  a check for validity of indices
nansignal(isnan(signal))=0;
nansignal(~isnan(signal))=1;


Cycle=[1:100]'*nan;
Cyclenan=[1:100]'*nan;
CycleLength=-100;
N=length(trigind)-1;
if N>1,
    for i=1:N,
        x=[trigind(i):trigind(i+1)]-trigind(i);
        CycleLength(i)=length(x);
        x=x*99/(trigind(i+1)-trigind(i));
        x=x+1;
        try
            Cycle(:,i)=interp1(x',signal(trigind(i):trigind(i+1))',[1:100]','cubic');
            Cyclenan(:,i)=interp1(x',nansignal(trigind(i):trigind(i+1))',[1:100]','cubic');
            
        catch
            Cyclenan(:,i)=nan;
            Cycle(:,i)=nan;
        end
    end
    Cyclenan(Cyclenan<1)=nan;
%     Cycle=Cycle.*Cyclenan;
    
    tmp=nanmean(Cycle(:,1:N)');
    Cycle(:,N+1)=tmp';
    tmp=nanstd(Cycle(:,1:N)');
    Cycle(:,N+2)=tmp';
    TimeGain=101/mean(CycleLength);
elseif N==1,
    x=[trigind(1):trigind(2)]-trigind(1);
    CycleLength=length(x);
    x=x.*99/(trigind(2)-trigind(1))+1;
    try
        Cycle(:,1)=interp1(x',signal(trigind(1):trigind(2))',[1:100]','cubic');
    catch
        Cycle(1:100,1)=nan;
    end
    TimeGain=100/CycleLength;
end


return
% ============================================
% END % ### NormalizeTimeBase