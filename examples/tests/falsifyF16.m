kf = modelF16();
% kf.spec = specification(halfspace([0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],0),'unsafeSet');
plot_vars = [12,13];
kf.verb=2;
kf.nResets=5;
kf.maxSims=1000;
kf.reach.split=true;

[kfSolns,allDatas] = falsify(kf);