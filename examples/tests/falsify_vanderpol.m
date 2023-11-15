kfModel = modelVanderpol();
plot_vars = [1,2];

[kfModel,trainset] = falsify(kfModel);

if kfModel.soln.falsified
    visualize_falsification(kfModel.soln.x, kfModel.soln.t, kfModel.spec, plot_vars)
    disp(['simulations required: ',num2str(kfModel.soln.sims)])
else
    disp("No falsifiying trace found")
end