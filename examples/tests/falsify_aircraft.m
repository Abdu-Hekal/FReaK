
% initialize seeds
rng(0)
pyrunfile("seed.py")
simsreq=[];
for j = 1:10
    kfModel = modelAircraftODE();

    [kfModel,trainset] = falsify(kfModel);

    if kfModel.soln.falsified
        disp(['simulations required: ',num2str(kfModel.soln.sims)])
        simsreq = [simsreq kfModel.soln.sims];
    else
        disp("No falsifiying trace found")
    end
end
avgSims=mean(simsreq);
medSims = median(simsreq);
fprintf('Avg Number of simulations: %.2f\n', avgSims);
fprintf('Median Number of simulations: %.2f\n', medSims);


