function repArch22()
% arch22ModelTransmission - runs all requirement formula for the
%  model transmission benchmark of the ARCH'22 falsification Category
%
% Syntax:
%   results = repArch22()
%
% Inputs:
%    -
%
% Outputs:
%    results -
%

% Author:       Abdelrahman Hekal
% Written:      23-Feb-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
model = model_AutoTransmission();

x = stl('x',3);
requirements = {; ...
    %     "AT1", globally(x(1) < 120,interval(0,20)); ...
    "AT2", globally(x(2) < 4750,interval(0,10)); ...
    %     "testAT2", globally(x(2) <= 4750,interval(0,10)); ...
%         "AT51", globally(implies(x(3)>=2 & finally(x(3)>=1 & x(3)<=1,interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    %      "AT6a", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
    %         "test", globally(x(1)<50 | x(1)>60,interval(10,30)),...
    %      "testAT6a", implies(globally(x(2)<3000,interval(0,4)),globally(x(1)<35,interval(0,4))); ...

    };

solns=dictionary(string.empty,cell.empty);
for i = 1:size(requirements, 1)
    for j = 1:1
        disp("--------------------------------------------------------")
        name = requirements{i, 1};
        eq = requirements{i, 2};

        model.spec = specification(eq,'logic');
        [model,~] = falsify(model);

        if j==1
            solns(name)={{model.soln}};
        else
            soln=solns(name);
            soln{1}{end+1}=model.soln;
            solns(name)=soln;
        end
    end
    avgKoopTime=getAvg(solns(name),'koopTime');
    avgMilpSetupTime=getAvg(solns(name),'milpSetupTime');
    avgMilpSolveTime=getAvg(solns(name),'milpSolvTime');
    avgRuntime=getAvg(solns(name),'runtime');
    avgTrain=getAvg(solns(name),'trainIter');
    avgFalsified=getAvg(solns(name),'falsified');
    %print info
    fprintf('Benchmark: %s\n', name);
    fprintf('Number of runs: %d\n', j);
    fprintf('Avg koopman time: %.2f seconds\n', avgKoopTime);
    fprintf('Avg milp setup time: %.2f seconds\n', avgMilpSetupTime);
    fprintf('Avg milp solve time: %.2f seconds\n', avgMilpSolveTime);
    fprintf('Avg total runtime: %.2f seconds\n', avgRuntime);
    fprintf('Avg training iterations: %.2f\n', avgTrain);
    [n,d] = numden(sym(avgFalsified));
    fprintf('Number of successful falsified traces: %d/%d\n', n,d);

end
end

function avg = getAvg(solns,metric)
    avg=0;
    for i=1:length(solns{1})
        avg=avg+solns{1}{i}.(metric);
    end
    avg=avg/i;
end

