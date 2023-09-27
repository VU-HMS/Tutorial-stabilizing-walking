function [vel,acc] = NumericalDerivative_CentralDiff(t, pos)
%NumericalDerivative Computes the first and second numerical derivative of
%the input vector pos (by fittings splines and computing the derivative of the spline). 
%   Input arguments:
%       (1): t = time vector
%       (2): pos = matrix with positions

% check if frames are on the first or second dim (pos)
[nr,nc] = size(pos);
if nc>nr
    BoolTranspose_Pos = true;
    pos = pos';
    nfr = nr;
else
    BoolTranspose_Pos = false;
    nfr = nc;
end
[nr,~] = size(t);
if nr ==1
    t = t';
end

% velocity
dy_dx = diff(pos) ./ diff(t);
vel = [dy_dx(1,:); 0.5 * (dy_dx(2:end,:) + dy_dx(1:end-1,:)); dy_dx(end,:)];
% acc
dy_dx = diff(vel) ./ diff(t);
acc = [dy_dx(1,:); 0.5 * (dy_dx(2:end,:) + dy_dx(1:end-1,:)); dy_dx(end,:)];


% transpose again for output if needed
if BoolTranspose_Pos
    vel = vel';
    acc = acc';
end






end