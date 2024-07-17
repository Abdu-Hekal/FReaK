function [Bdata,phi,rob]=bReachRob(coraSpec,t,x,u)
%define variables using state variables
vars=arrayfun(@(i) ['x', num2str(i)], 1:size(x,2), 'UniformOutput', false)';
if ~isempty(u)
    trace = [t,x,u'];
    % add inputs to variables
    uvars=arrayfun(@(i) ['u', num2str(i)], 1:size(u,1), 'UniformOutput', false)';
    vars=[vars;uvars];
else
    trace = [t,x];
end

Bdata = BreachTraceSystem(vars');
Bdata.AddTrace(trace);
stl=coraBreachConvert(coraSpec.set);
phi = STL_Formula('phi',stl);

% Rphi = BreachRequirement(phi);
% rob=Rphi.Eval(Bdata);

%faster robustness evaluation than Rphi.Eval
traj=Bdata.P.traj{1};
val = STL_Eval(Bdata.Sys, phi, Sselect(Bdata.P,1), traj, traj.time);
rob=val(1);


end