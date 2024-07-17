function enhancement()
% enhancement - runs experiments in paper on enhancements
%
% Syntax:
%   results = enhancement()
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

bench.kfModel = @modelNeural;
x = stl('x',2);
u = stl('u',1);
alpha=0.005;
beta=0.03;
beta2=0.04;
bench.requirements = {; ...
    "NN", 40/12, globally(implies(abs(x(1)-u(1))>alpha+beta*abs(u(1)),finally(globally(~(alpha+beta*abs(u(1))<=abs(x(1)-u(1))),interval(0,1)),interval(0,2))),interval(1,37)); ...
    "NN2", 40/12, globally(implies(abs(x(1)-u(1))>alpha+beta2*abs(u(1)),finally(globally(~(alpha+beta2*abs(u(1))<=abs(x(1)-u(1))),interval(0,1)),interval(0,2))),interval(1,37)); ...
    "NNx", 40/12, (finally(x(1)>3.2,interval(0,1))) & (finally(globally(x(1)>1.75 & x(1)<2.25,interval(0,0.5)),interval(1,1.5))) & (globally(x(1)>1.825 & x(1)<2.175,interval(2,3))) ; ...
    };
benches{end+1} = bench;

bench.kfModel = @modelSC;
x = stl('x',4);
bench.requirements = {; ...
    "SC",0.1, globally(x(4)>=87 & x(4)<=87.5,interval(30,35)) ; ...
    };
benches{end+1} = bench;

% Start recording the command line output to a file
diary('enhancement.txt');

enhancements=[0,1,2,3,4];
for e = 1:numel(enhancements)
    solns=dictionary(string.empty,cell.empty);
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
            kfModel.timeout=1000;
            if enhancements(e)==0
                kfModel.nResets=inf; kfModel.rmRand=0; kfModel.offsetStrat=0; kfModel.reach.on=0;
            elseif enhancements(e)==1
                kfModel.nResets=5; kfModel.rmRand=0; kfModel.offsetStrat=0; kfModel.reach.on=0;
            elseif enhancements(e)==2
                kfModel.nResets=5; kfModel.rmRand=0; kfModel.offsetStrat=-1; kfModel.reach.on=0;
            elseif enhancements(e)==3
                kfModel.nResets=5; kfModel.rmRand=0; kfModel.offsetStrat=-1; kfModel.reach.on=1;
            elseif enhancements(e)==4
                kfModel.nResets=5; kfModel.rmRand=1; kfModel.offsetStrat=-1; kfModel.reach.on=1;
            end

            %skip at6 benchmarks if no offset (unable to falsify any of them without offset)
            if enhancements(e)<2 && contains(name,'AT6')
                break;
            end

            if name == "NNx"
                kfModel.U = interval(1.95,2.05);
            end

            kfModel.spec = specification(eq,'logic');
            kfModel.runs=10;
            kfSolns = falsify(kfModel);

            %print info
            fprintf('Benchmark: %s\n', name);
            fprintf('enhancement=%d \n',enhancements(e));
            printInfo(kfSolns)


        end
    end
end
% Stop recording the command line output
diary off;
end

