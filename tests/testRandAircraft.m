% initialize seeds
rng(0)
pyrunfile("seed.py")
kfModel = modelAircraftODE();

minRob=inf;
tsim = (0:kfModel.dt:kfModel.T)'; %define time points for interpolating simulation
for i=1:100
    [t, x, u] = sampleSimulation(kfModel);
    usim = interp1(u(:,1),u(:,2:end),tsim,kfModel.inputInterpolation,"extrap"); %interpolate and extrapolate input points
    usim =  max(kfModel.U.inf',min(kfModel.U.sup',usim)); %ensure that extrapolation is within input bounds
    usim = [tsim,usim];
    interpX = interp1(t,x,tsim,kfModel.trajInterpolation); %interpolate trajectory at granulated time points for checking correctness

     spec=kfModel.spec(1);
     [~,~,robustness] = bReachRob(spec,tsim,interpX,usim(:,2:end)');
     if robustness < minRob
         minRob=robustness;
     end
     if robustness<0
         break;
     end
end
i
minRob