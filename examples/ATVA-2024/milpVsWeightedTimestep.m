function milpVsWeightedTimestep()
% timestep - runs experiments in paper on timestep
%
% Syntax:
%   results = timestep()
%
% Inputs:
%    -
%
% Outputs:
%    results -
%

% Author:       Abdelrahman Hekal
% Written:      30-Oct-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
benches = {}; %empty cell to store benchmarks
% Chasing cars benchmark
bench.kfModel = @modelCars;
x = stl('x',5);
bench.requirements = {; ...
    "CC5", 100/50, globally(finally(implies(globally(x(2)-x(1)>=9,interval(0,5)),globally(x(5)-x(4)>=9,interval(5,20))),interval(0,8)),interval(0,72)); 
    };
benches{end+1} = bench;

% Start recording the command line output to a file
diary('fig4.txt');

%50,40,30,20,10 control points
timesteps=[100/10,100/20,100/30,100/40,100/50];
%weighted approach true or false
weighted=[0,1];
for w=1:numel(weighted)
    for t = 1:numel(timesteps)
        for b = 1:length(benches)
            bench = benches{b};
            req = bench.requirements;
            for i = 1:size(req, 1)
                % initialize seeds
                rng(0)
                pyrunfile("seed.py")
                disp("--------------------------------------------------------")
                kfModel = bench.kfModel();
                name = req{i, 1};
                eq = req{i, 3};
                %settings
                kfModel.ak.dt=timesteps(t);

                kfModel.spec = specification(eq,'logic');
                kfModel.runs=10;
                kfModel.verb=0;
                kfModel.nResets=5; %'auto'
                kfModel.resetStrat=0;
                kfModel.maxSims=1000;

                %turn on.off weighted approach, if off: MILP is used 
                kfModel.solver.autoAddTimePoints=true*weighted(w);
                kfModel.solver.autoAddConstraints=2*weighted(w);

                kfSolns = falsify(kfModel);

                %print info
                fprintf('\nBenchmark: %s\n', name);
                fprintf('Timestep=%.1f \n',timesteps(t));
                printInfo(kfSolns)
            end
        end
    end
end
% Stop recording the command line output
diary off;
end

