function [Bdata,phi,rob]=bReachRob(coraSpec,t,x,u)
%define variables using trajectory and inputs
xvars=arrayfun(@(i) ['x', num2str(i)], 1:size(x,2), 'UniformOutput', false)';
uvars=arrayfun(@(i) ['u', num2str(i)], 1:size(u,1), 'UniformOutput', false)';
vars=[xvars;uvars];

Bdata = BreachTraceSystem(vars');
trace = [t,x,u'];
Bdata.AddTrace(trace);
stl=replace(coraBlustlConvert(coraSpec.set),"(t)","[t]");
phi = STL_Formula('phi',stl);
Rphi = BreachRequirement(phi);

rob=Rphi.Eval(Bdata);
end