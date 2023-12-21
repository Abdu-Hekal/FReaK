%set seed
% rng(0)
% pyrunfile("seed.py")

kf = modelF16();
kf.verb=2;

kf.maxSims=100;
kf.nResets=10;
bestSoln.rob=inf;

%partition non-exact dimensions
ind=find(kf.R0.sup-kf.R0.inf);
R0_ = partition(kf.R0(ind),4);

for ii=1:numel(R0_)
    kf.R0(ind)=R0_{ii};
    soln=falsify(kf);
    if soln.best.rob < bestSoln.rob
        bestSoln.rob = soln.best.rob;
        bestSoln.R0=kf.R0;
    end
end

fprintf('Best Robustness=%.3f\n',bestSoln.rob)

