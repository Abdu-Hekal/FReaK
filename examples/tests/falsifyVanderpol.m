function falsifyVanderpol()

kfModel = modelVanderpol();
plot_vars = [1,2];

soln = falsify(kfModel);

if soln.falsified
    visualizeFalsification(soln.best.x, soln.best.t, kfModel.spec, plot_vars)
    disp(['simulations required: ',num2str(soln.sims)])
else
    disp("No falsifiying trace found")
end

end