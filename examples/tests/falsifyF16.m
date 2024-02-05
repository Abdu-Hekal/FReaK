kf = modelF16();
% kf.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];
kf.verb=2;
kf.nResets=10;
kf.maxSims=1;

[kfSolns,allDatas] = falsify(kf);
kfSoln=kfSolns{1};
allData=allDatas{1};
visualizeTrain(allData,kf.ak.dt, plot_vars)