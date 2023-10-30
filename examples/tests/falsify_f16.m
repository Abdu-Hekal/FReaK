% load your python environment that autokoopman is installed & imported in
% pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
% py.importlib.import_module('autokoopman');

kfModel = modelF16();
% kfModel.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];

[kfModel,trainset] = falsify(kfModel);

if kfModel.soln.falsified
    visualize_falsification(kfModel.soln.x, trainset.t{1}, kfModel.spec, plot_vars)
    disp(['simulations required: ',num2str(kfModel.soln.sims)])
    visualize_train(trainset, plot_vars)
else
    disp("No falsifiying trace found")
end

visualize_train(trainset, plot_vars)