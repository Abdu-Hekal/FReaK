%set seed
rng(12)

x = stl('x',3);
kfModel = modelAutoTransmission();

eq=globally(implies(~(x(3)>=1 & x(3)<=1) & finally(x(3)>=1 & x(3)<=1,interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30));
conjunctiveeq = conjunctiveNormalForm(eq);

spec = specification(eq,'logic');
conjunctiveSpec=specification(conjunctiveeq,'logic');

[tout, yout, u]=kfModel.sampleSimulation();
tsim = (0:kfModel.dt:kfModel.T)'; %define time points for interpolating simulation
usim = interp1(u(:,1),u(:,2:end),tsim,kfModel.inputInterpolation,"extrap"); %interpolate and extrapolate input points
usim =  max(kfModel.U.inf',min(kfModel.U.sup',usim)); %ensure that extrapolation is within input bounds
usim = [tsim,usim];
interpX = interp1(tout,yout,tsim,kfModel.trajInterpolation); %interpolate trajectory at granulated time points for checking correctness


save('tsim.mat','tsim')
save('usim.mat','usim')
save('interpX.mat','interpX')


[~,~,rob] = bReachRob(spec,tsim,interpX,usim(:,2:end)');
[~,~,conjunctRob] = bReachRob(conjunctiveSpec,tsim,interpX,usim(:,2:end)');


rob
conjunctRob

