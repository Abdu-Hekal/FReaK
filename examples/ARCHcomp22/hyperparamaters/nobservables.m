function nobservables()
% nobservables - runs experiments in paper on number of observables
%
% Syntax:
%   results = nobservables()
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
%Model transmission benchmark
bench.kfModel = @modelAutoTransmission;
x = stl('x',3);
bench.requirements = {; ...
    "AT1", 1, globally(x(1) < 120,interval(0,20)); ...
    };
benches{end+1} = bench;
% Chasing cars benchmark
bench.kfModel = @modelCars;
x = stl('x',5);
bench.requirements = {; ...
    "CC1", 10, globally(x(5)-x(4)<=40,interval(0,100)); ...
    };
benches{end+1} = bench;

% Start recording the command line output to a file
diary('nobservables.txt');

observables=[10,20,30,40,50];
for o = 1:numel(observables)
    for reach=0:1
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
                kfModel.ak.dt = req{i, 2};
                eq = req{i, 3};
                %settings
                kfModel.reach.on=reach;
                kfModel.ak.nObs=observables(o);
                kfModel.ak.rank=[1,observables(o),observables(o)/5];

                if name == "NNx"
                    kfModel.U = interval(1.95,2.05);
                end
                if name == "AFC33"
                    kfModel.U = interval([61.2;900],[81.2;1100]);
                end

                kfModel.spec = specification(eq,'logic');
                kfModel.runs=10;
                kfSolns = falsify(kfModel);

                %print info
                fprintf('Benchmark: %s\n', name);
                if reach==0
                    disp('- Baseline -')
                else
                    disp('- Reach -')
                end
                fprintf('Number of observables=%d \n',observables(o));
                printInfo(kfSolns)
            end
        end
    end
end
% Stop recording the command line output
diary off;
end


