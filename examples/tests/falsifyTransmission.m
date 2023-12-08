
% initialize seeds
% rng(18)
% pyrunfile("seed.py")

kfModel = modelAutoTransmission();
plot_vars = [1,2];

x = stl('x',3);
eq = implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4)));
kfModel.spec = specification(eq,'logic');

kfModel.maxSims=100;
kfModel.verb=3;
kfModel.nResets='auto';

[kfSoln,trainset] = falsify(kfModel);

if kfSoln.falsified
    disp(['simulations required: ',num2str(kfSoln.sims)])
else
    disp("No falsifiying trace found")
end

visualizeTrain(trainset, kfSoln.koopModel, plot_vars,'Speed','Angular velocity')
%settings for figure
% figure_settings(gcf);
% export_fig training.pdf