function repArch22()
% arch22ModelTransmission - runs all requirement formula for the
%  benchmarks of the ARCH'22 falsification Category
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
benches = {}; %empty cell to store benchmarks
%Model transmission benchmark
% bench.kfModel = modelAutoTransmission();
% x = stl('x',3);
% bench.requirements = {; ...
% %             "AT1", globally(x(1) < 120,interval(0,20)); ...
% %         "AT2", globally(x(2) < 4750,interval(0,10)); ...
% %             "AT51", globally(implies(~(x(3)>=1 & x(3)<=1) & finally(x(3)>=1 & x(3)<=1,interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
% %             "AT52", globally(implies(~(x(3)>=2 & x(3)<=2) & finally(x(3)>=2 & x(3)<=2,interval(0.001,0.1)),finally(globally(x(3)>=2 & x(3)<=2,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
% %             "AT53", globally(implies(~(x(3)>=3 & x(3)<=3) & finally(x(3)>=3 & x(3)<=3,interval(0.001,0.1)),finally(globally(x(3)>=3 & x(3)<=3,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
% %             "AT54", globally(implies(~(x(3)>=4 & x(3)<=4) & finally(x(3)>=4 & x(3)<=4,interval(0.001,0.1)),finally(globally(x(3)>=4 & x(3)<=4,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
% 
% %     "AT6a", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
%     %          "AT6b", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))); ...
%     %          "AT6c", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20))); ...
% 
% %             "test", globally(x(1)<50 | x(1)>60,interval(10,30)),...
%          "test", finally(x(2)<3000 | x(2)>4000,interval(10,30)),...
%     };
% benches{end+1} = bench;
% Chasing cars benchmark
bench.kfModel = modelCars();
bench.kfModel.trainRand=0;
x = stl('x',5);
bench.requirements = {; ...
    "CC1", globally(x(5)-x(4)<=40,interval(0,100)); ...
    };
benches{end+1} = bench;

solns=dictionary(string.empty,cell.empty);
for b = 1:length(benches)
    bench = benches{1};
    req = bench.requirements;
    kfModel = bench.kfModel;
    for i = 1:size(req, 1)
        for j = 1:10
            disp("--------------------------------------------------------")
            name = req{i, 1};
            eq = req{i, 2};

            kfModel.spec = specification(eq,'logic');
            [kfModel,~] = falsify(kfModel);

            if j==1
                solns(name)={{kfModel.soln}};
            else
                soln=solns(name);
                soln{1}{end+1}=kfModel.soln;
                solns(name)=soln;
            end
            fprintf('falsified iteration %d \n',j);
        end
        avgKoopTime=mean(getMetrics(solns(name),'koopTime'));
        avgMilpSetupTime=mean(getMetrics(solns(name),'milpSetupTime'));
        avgMilpSolveTime=mean(getMetrics(solns(name),'milpSolvTime'));
        avgRuntime=mean(getMetrics(solns(name),'runtime'));
        sims=getMetrics(solns(name),'sims');
        avgSims=mean(sims);
        medianSims=median(sims);
        avgFalsified=sum(getMetrics(solns(name),'falsified'));
        %print info
        fprintf('Benchmark: %s\n', name);
        fprintf('Number of runs: %d\n', j);
        fprintf('Avg koopman time: %.2f seconds\n', avgKoopTime);
        fprintf('Avg milp setup time: %.2f seconds\n', avgMilpSetupTime);
        fprintf('Avg milp solve time: %.2f seconds\n', avgMilpSolveTime);
        fprintf('Avg total runtime: %.2f seconds\n', avgRuntime);
        fprintf('Avg Number of simulations: %.2f\n', avgSims);
        fprintf('Median Number of simulations: %.2f\n', medianSims);
        fprintf('Number of successful falsified traces: %d/%d\n', avgFalsified,j);

    end
    save("solns.mat","solns")
end
end

function list = getMetrics(solns,metric)
list=[];
for i=1:length(solns{1})
    list=[list,solns{1}{i}.(metric)];
end
end
