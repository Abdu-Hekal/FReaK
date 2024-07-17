function timestep()
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
%Model transmission benchmark
bench.kfModel = @modelAutoTransmission;
x = stl('x',3);
bench.requirements = {; ...
    "AT1", 1, globally(x(1) < 120,interval(0,20)); ...
    "AT2", 1,globally(x(2) < 4750,interval(0,10)); ...
    "AT51", 1,globally(implies(~(x(3)>=1 & x(3)<=1) & finally(x(3)>=1 & x(3)<=1,interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT52", 1, globally(implies(~(x(3)>=2 & x(3)<=2) & finally(x(3)>=2 & x(3)<=2,interval(0.001,0.1)),finally(globally(x(3)>=2 & x(3)<=2,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT53", 1, globally(implies(~(x(3)>=3 & x(3)<=3) & finally(x(3)>=3 & x(3)<=3,interval(0.001,0.1)),finally(globally(x(3)>=3 & x(3)<=3,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT54", 1, globally(implies(~(x(3)>=4 & x(3)<=4) & finally(x(3)>=4 & x(3)<=4,interval(0.001,0.1)),finally(globally(x(3)>=4 & x(3)<=4,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT6a", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
    "AT6b", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))); ...
    "AT6c", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20))); ...
    "AT6abc", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20)))
    };
benches{end+1} = bench;
% Chasing cars benchmark
bench.kfModel = @modelCars;
x = stl('x',5);
bench.requirements = {; ...
    "CC1", 10, globally(x(5)-x(4)<=40,interval(0,100)); ...
    "CC2", 10, globally(finally(x(5)-x(4)>=15,interval(0,30)),interval(0,70)); ...
    "CC3", 10, globally(globally(x(2)-x(1)<=20,interval(0,20)) | finally(x(5)-x(4)>=40,interval(0,20)),interval(0,80)); ...
    "CC4", 10, globally(finally(globally(x(5)-x(4)>=8,interval(0,20)),interval(0,30)),interval(0,65)); ...
    "CC5", 10, globally(finally(implies(globally(x(2)-x(1)>=9,interval(0,5)),globally(x(5)-x(4)>=9,interval(5,20))),interval(0,8)),interval(0,72)); ...
    "CCx", 10, globally(x(2)-x(1)>7.5,interval(0,50)) & globally(x(3)-x(2)>7.5,interval(0,50)) & globally(x(4)-x(3)>7.5,interval(0,50)) & globally(x(5)-x(4)>7.5,interval(0,50)); ...
    };
benches{end+1} = bench;

bench.kfModel = @modelSC;
x = stl('x',4);
bench.requirements = {; ...
    %     "SC",0.1, globally(x(4)>=87 & x(4)<=87.5,interval(30,35)) ; ...
    };
benches{end+1} = bench;

% Start recording the command line output to a file
diary('timestep.txt');

timesteps=[0.1,0.5,1,2.5,5,10];
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
            kfModel.timeout=1000;
            kfModel.ak.dt=timesteps(t);

            if name == "NNx"
                kfModel.U = interval(1.95,2.05);
            end

            kfModel.spec = specification(eq,'logic');
            kfModel.runs=10;
            kfSolns = falsify(kfModel);

            %print info
            fprintf('Benchmark: %s\n', name);
            fprintf('Timestep=%.1f \n',timesteps(t));
            printInfo(kfSolns)
        end
    end
end
% Stop recording the command line output
diary off;
end

