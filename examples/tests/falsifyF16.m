% load your python environment that autokoopman is installed & imported in
% pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
% py.importlib.import_module('autokoopman');

kf = modelF16();
% kf.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];
kf.verb=2;
kf.nResets=10;

[kf,trainset] = falsify(kf);

if kf.soln.falsified
    visualizeFalsification(kf.soln.x, trainset.t{1}, kf.spec, plot_vars)
    disp(['simulations required: ',num2str(kf.soln.sims)])
    visualize_train(trainset, plot_vars)
else
    disp("No falsifiying trace found")
end

visualize_train(trainset, plot_vars)