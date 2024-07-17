%set seed
% rng(0)
% pyrunfile("seed.py")

kf = modelF16();
kf.verb=2;

kf.maxSims=100;
kf.nResets=10;
bestSoln.rob=inf;
dim=6;

bestSoln=recursiveSplit(kf,bestSoln,dim,1,9);
fprintf('Best Robustness=%.2f, on level %d \n',bestSoln.rob,bestSoln.level)
fprintf('Best interval: [%.3f, %.3f] \n',bestSoln.R0.inf(dim), bestSoln.R0.sup(dim))


function bestSoln=recursiveSplit(kf,bestSoln,dim, currentLevel, maxLevel)
    if currentLevel > maxLevel
%         disp(['Reached maximum level: ', num2str(maxLevel)]);
        return;
    end

    %falsify and update best soln
    fprintf('split level %d \n',currentLevel)
    disp('--->')
    solns = falsify(kf);
    soln=solns{1};
    if soln.best.rob < bestSoln.rob
        bestSoln.rob = soln.best.rob;
        bestSoln.level=currentLevel;
        bestSoln.R0=kf.R0;
    end

    % Split the initial set into two halves and falsify each half
    R0_ = kf.R0.split(dim);
    for sp=1:numel(R0_)
        kf.R0 = R0_{sp};
        bestSoln=recursiveSplit(kf,bestSoln,dim, currentLevel + 1, maxLevel);
    end
end
