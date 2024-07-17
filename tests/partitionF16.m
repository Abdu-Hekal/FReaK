%set seed
% rng(0)
% pyrunfile("seed.py")

kf = modelF16();
kf.verb=1;

kf.maxSims=100;
kf.nResets=10;
bestSoln.rob=inf;

%partition non-exact dimensions
ind=find(kf.R0.sup-kf.R0.inf);
R0_ = partition(kf.R0(ind),10);

for ii=1:numel(R0_)
    kf.R0(ind)=R0_{ii};
    solns=falsify(kf);
    soln=solns{1};
    fprintf('Best Robustness=%.3f\n',soln.best.rob)
    if soln.best.rob < bestSoln.rob
        disp('New best Robustness found!')
        bestSoln.rob = soln.best.rob;
        bestSoln.R0=kf.R0;
    end
    if bestSoln.rob < 0
        break;
    end
end

fprintf('Best Robustness=%.3f\n',bestSoln.rob)
save('bestSoln.mat',"bestSoln")

