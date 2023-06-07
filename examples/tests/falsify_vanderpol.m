% load your python environment that autokoopman is installed & imported in
% pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
% py.importlib.import_module('autokoopman');

kfModel = modelVanderpol();
kfModel.maxTrainSize=10; %maximum number of training trajectories before quitting
plot_vars = [1,2];

[kfModel,trainset] = falsify(kfModel);

if kfModel.soln.falsified
    visualize_falsification(kfModel.soln.x, trainset.t{1}, kfModel.spec, plot_vars)
    disp(['training iterations required: ',num2str(kfModel.soln.trainIter)])
    visualize_train(trainset, plot_vars)
else
    disp("No falsifiying trace found")
end

visualize_train(trainset, plot_vars)