function falsifyVanderpol()

kfModel = modelVanderpol();
kfModel.verb=3;
plot_vars = [1,2];

solns = falsify(kfModel);
soln=solns{1};

if soln.falsified
    visualizeFalsification(soln.best.x, soln.best.t, kfModel.spec, plot_vars)
end

end