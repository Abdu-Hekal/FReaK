% load your python environment that autokoopman is installed & imported in
% pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
% py.importlib.import_module('autokoopman');

model = model_vanderpol();
model.max_train_size=10; %maximum number of training trajectories before quitting
plot_vars = [1,2];

[model,trainset] = falsify(model);

if falsified
    disp("falsifying trace found")
    disp(['training iterations required: ',num2str(train_iter)])
    visualize_falsification(model.soln.x, trainset.t{1}, model.spec, plot_vars)
else
    disp("No falsifiying trace found")
end

visualize_train(trainset, plot_vars)