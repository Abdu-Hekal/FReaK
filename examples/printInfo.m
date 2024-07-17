function printInfo(solns)
% soln: cells of kfModel.soln, see repArch22 for example

avgFalsified=sum(getMetrics(solns,'falsified'));
fprintf('Number of successful falsified traces: %d/%d\n', avgFalsified,numel(solns));

if avgFalsified>0
    avgKoopTime=mean(getMetrics(solns,'koopTime'));
    avgOptimizerTime=mean(getMetrics(solns,'optimTime'));
    avgSimTime=mean(getMetrics(solns,'simTime'));
    avgRuntime=mean(getMetrics(solns,'runtime'));
    sims=getMetrics(solns,'sims');
    avgSims=mean(sims);
    medianSims=median(sims);
    upper_quart = upper_quartile(sims);
    lower_quart = lower_quartile(sims);
    upper_whisk = upper_whisker(sims);
    lower_whisk = lower_whisker(sims);
    outliers = comp_outliers(sims, upper_whisk, lower_whisk);

    disp('Computation Time -->');
    fprintf('Avg koopman time: %.2f seconds\n', avgKoopTime);
    fprintf('Avg optimization time: %.2f seconds\n', avgOptimizerTime);
    fprintf('Avg simulation time: %.2f seconds\n', avgSimTime);
    fprintf('Avg total runtime: %.2f seconds\n', avgRuntime);
    fprintf('R: %.2f\n',100*avgSimTime/avgRuntime)
    disp('Simulations -->');
    fprintf('average: %.2f\n', avgSims);
    fprintf('median: %.2f\n', medianSims);
    fprintf('upper quartile= %.2f\n', upper_quart);
    fprintf('lower quartile= %.2f\n', lower_quart);
    fprintf('upper whisker= %.2f\n', upper_whisk);
    fprintf('lower whisker= %.2f\n', lower_whisk);
    fprintf('Outliers=[')
    for i=1:length(outliers)
        fprintf(' %d,', outliers(i));
    end
    fprintf(']\n');

end

end

function list = getMetrics(solns,metric)
list=[];
for i=1:length(solns)
    if solns{i}.falsified
        list=[list,solns{i}.(metric)];
    end
end
end

function res=upper_quartile(data)
res= prctile(data, 75);
end

function res=lower_quartile(data)
res= prctile(data, 25);
end

function res=upper_whisker(data)
iqr = upper_quartile(data) - lower_quartile(data);
res= min(upper_quartile(data) + 1.5 * iqr, max(data));
end

function res=lower_whisker(data)
iqr = upper_quartile(data) - lower_quartile(data);
res= max(lower_quartile(data) - 1.5 * iqr, min(data));
end

function res=comp_outliers(data, upper_whisk, lower_whisk)
res = [];
for ii=1:numel(data)
    if data(ii) > upper_whisk || data(ii) < lower_whisk
        res = [res data(ii)];
    end
end
end
