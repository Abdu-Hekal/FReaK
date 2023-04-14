function rob=bReachRob(coraSpec,x,t)
vars=coraSpec.set.getVariables;
Bdata = BreachTraceSystem(vars');
trace = [t,x];
Bdata.AddTrace(trace);
stl=replace(coraBlustlConvert(coraSpec.set),"(t)","[t]");
phi = STL_Formula('phi',stl);
Rphi = BreachRequirement(phi);


robustness = @(Bdata) Rphi.Eval(Bdata);

rob=robustness(Bdata);
end