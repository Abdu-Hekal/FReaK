function [Bdata,phi,rob]=bReachRob(coraSpec,x,t)
vars=coraSpec.set.getVariables;
Bdata = BreachTraceSystem(vars');
trace = [t,x];
Bdata.AddTrace(trace);
stl=replace(coraBlustlConvert(coraSpec.set),"(t)","[t]");
phi = STL_Formula('phi',stl);
Rphi = BreachRequirement(phi);

rob=Rphi.Eval(Bdata);
end