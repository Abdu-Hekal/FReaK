% load your python environment that autokoopman is installed & imported in
% pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
% py.importlib.import_module('autokoopman');

model = model_AutoTransmission();
plot_vars = [1,2];

[model,trainset] = falsify(model);

if model.soln.falsified
    disp(" ")
    disp("falsifying trace found!")
    visualize_falsification(model.soln.x, trainset.t{1}, model.spec, plot_vars)
else
    disp("No falsifying trace found!")
end
disp(['training iterations: ',num2str(train_iter)])
visualize_train(trainset, plot_vars)