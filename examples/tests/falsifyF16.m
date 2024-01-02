kf = modelF16();
% kf.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];
kf.verb=2;
kf.nResets=10;
kf.maxSims=1;

[kfSolns,trainset] = falsify(kf);
kfSoln=kfSolns{1};
if kfSoln.falsified
    visualizeFalsification(kfSoln.best.x, trainset.t{1}, kf.spec, plot_vars)
end

visualizeTrain(trainset,kfSoln.koopModel, plot_vars)