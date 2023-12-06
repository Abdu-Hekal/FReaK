kf = modelF16();
% kf.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];
kf.verb=2;
kf.nResets=10;
kf.maxSims=1;

[kfSoln,trainset] = falsify(kf);

if kfSoln.falsified
    visualizeFalsification(kfSoln.best.x, trainset.t{1}, kf.spec, plot_vars)
    disp(['simulations required: ',num2str(kf.soln.sims)])
else
    disp("No falsifiying trace found")
end

visualizeTrain(trainset,kfSoln.koopModel, plot_vars)