function [offsetCount] = getRobOffset(obj,x,t,rob)
% getRobOffset - compute which subformula to offset in stl to falsify
% system
%
% Syntax:
%    rho = getRobOffset(obj,x,t,rob)
%
% Inputs:
%    obj - logic formula (class stl)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%    rob - current robustness of trace 
%
% Outputs:
%    offsetCount - count corresponding to stl subformula to offset
offsetCount=0;
offset=rob;
while rob > 0
    rob = analyseRobustness(obj,x,t,offset,offsetCount);
    offsetCount=offsetCount+1;
end
end