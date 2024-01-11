
% initialize seeds
% rng(18)
% pyrunfile("seed.py")

kfModel = modelAutoTransmission();
plot_vars = [1,2];

x = stl('x',3);
eq = implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4)));
kfModel.spec = specification(eq,'logic');

kfModel.maxSims=1;
kfModel.verb=2;

[kfSolns,allData] = falsify(kfModel);
kfSoln=kfSolns{1};

if kfSoln.falsified
    disp(['simulations required: ',num2str(kfSoln.sims)])
else
    disp("No falsifiying trace found")
end

visualizeTrain(allData,kfModel.ak.dt,plot_vars,'Speed','Angular velocity')
%settings for figure
% figure_settings(gcf);
% export_fig training.pdf