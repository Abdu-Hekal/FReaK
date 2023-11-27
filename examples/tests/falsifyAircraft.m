
% initialize seeds
rng(0)
pyrunfile("seed.py")
simsreq=[];

solns={};
for j = 1:10
    kfModel = modelAircraftODE();

    [kfModel,trainset] = falsify(kfModel);

    if kfModel.soln.falsified
        disp(['simulations required: ',num2str(kfModel.soln.sims)])
        simsreq = [simsreq kfModel.soln.sims];
        solns{end+1}=kfModel.soln;
    else
        disp("No falsifiying trace found")
    end
end
avgKoopTime=mean(getMetrics(solns,'koopTime'));
avgReachTime=mean(getMetrics(solns,'reachTime'));
avgMilpSetupTime=mean(getMetrics(solns,'milpSetupTime'));
avgMilpSolveTime=mean(getMetrics(solns,'milpSolvTime'));
avgSimTime=mean(getMetrics(solns,'simTime'));
avgRuntime=mean(getMetrics(solns,'runtime'));
sims=getMetrics(solns,'sims');
avgSims=mean(sims);
medianSims=median(sims);
avgFalsified=sum(getMetrics(solns,'falsified'));
%print info
fprintf('Benchmark: Aircraft');
fprintf('Number of runs: %d\n', j);
fprintf('Avg koopman time: %.2f seconds\n', avgKoopTime);
fprintf('Avg Reachability time: %.2f seconds\n', avgReachTime);
fprintf('Avg milp setup time: %.2f seconds\n', avgMilpSetupTime);
fprintf('Avg milp solve time: %.2f seconds\n', avgMilpSolveTime);
fprintf('Avg simulation time: %.2f seconds\n', avgSimTime);
fprintf('Avg total runtime: %.2f seconds\n', avgRuntime);
fprintf('Number of successful falsified traces: %d/%d\n', avgFalsified,j);
fprintf('Avg Number of simulations: %.2f\n', avgSims);
fprintf('Median Number of simulations: %.2f\n', medianSims);
fprintf('R: %.2f\n',100*avgSimTime/avgRuntime)


function list = getMetrics(solns,metric)
list=[];
for i=1:length(solns)
    list=[list,solns{i}.(metric)];
end
end



