function repArch24_2()
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
bench.kfModel = @modelAutoTransmission2;
x = stl('x',3);
bench.requirements = {; ...
    "AT","AT1", 5, globally(x(1) < 120,interval(0,20)); ...
    "AT","AT2", 5,globally(x(2) < 4750,interval(0,10)); ...
    "AT","AT51", 5,globally(implies(~(x(3)>=1 & x(3)<=1) & finally(x(3)>=1 & x(3)<=1,interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT","AT52", 5, globally(implies(~(x(3)>=2 & x(3)<=2) & finally(x(3)>=2 & x(3)<=2,interval(0.001,0.1)),finally(globally(x(3)>=2 & x(3)<=2,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT","AT53", 5, globally(implies(~(x(3)>=3 & x(3)<=3) & finally(x(3)>=3 & x(3)<=3,interval(0.001,0.1)),finally(globally(x(3)>=3 & x(3)<=3,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT","AT54", 5, globally(implies(~(x(3)>=4 & x(3)<=4) & finally(x(3)>=4 & x(3)<=4,interval(0.001,0.1)),finally(globally(x(3)>=4 & x(3)<=4,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    "AT","AT6a", 5, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
    "AT","AT6b", 5, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))); ...
    "AT","AT6c", 5, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20))); ...
    "AT","AT6abc", 5, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20)))
    };
benches{end+1} = bench;
% Chasing cars benchmark
bench.kfModel = @modelCars2;
x = stl('x',5);
bench.requirements = {; ...
    "CC","CC1", 5, globally(x(5)-x(4)<=40,interval(0,100)); ...
    "CC","CC2", 5, globally(finally(x(5)-x(4)>=15,interval(0,30)),interval(0,70)); ...
    "CC","CC3", 5, globally(globally(x(2)-x(1)<=20,interval(0,20)) | finally(x(5)-x(4)>=40,interval(0,20)),interval(0,80)); ...
    "CC","CC4", 5, globally(finally(globally(x(5)-x(4)>=8,interval(0,5)),interval(0,30)),interval(0,65)); ...
    "CC","CC5", 5, globally(finally(implies(globally(x(2)-x(1)>=9,interval(0,5)),globally(x(5)-x(4)>=9,interval(5,20))),interval(0,8)),interval(0,72)); ...
    "CC","CCx", 5, globally(x(2)-x(1)>7.5,interval(0,50)) & globally(x(3)-x(2)>7.5,interval(0,50)) & globally(x(4)-x(3)>7.5,interval(0,50)) & globally(x(5)-x(4)>7.5,interval(0,50)); ...
    };
benches{end+1} = bench;

bench.kfModel = @modelNeural2;
x = stl('x',2);
u = stl('u',1);
alpha=0.005;
beta=0.03;
beta2=0.04;
bench.requirements = {; ...
    "NN","(NN 0.005 0.03)", 40/3, globally(implies(abs(x(1)-u(1))>alpha+beta*abs(u(1)),finally(globally(~(alpha+beta*abs(u(1))<=abs(x(1)-u(1))),interval(0,1)),interval(0,2))),interval(1,37)); ...
    "NN","(NN 0.005 0.04)", 40/3, globally(implies(abs(x(1)-u(1))>alpha+beta2*abs(u(1)),finally(globally(~(alpha+beta2*abs(u(1))<=abs(x(1)-u(1))),interval(0,1)),interval(0,2))),interval(1,37)); ...
    "NN","NNx", 1, (finally(x(1)>3.2,interval(0,1))) & (finally(globally(x(1)>1.75 & x(1)<2.25,interval(0,0.5)),interval(1,1.5))) & (globally(x(1)>1.825 & x(1)<2.175,interval(2,3))) ; ...
    };
benches{end+1} = bench;

bench.kfModel = @modelSC2;
x = stl('x',4);
bench.requirements = {; ...
    %     "SC",0.05, globally(x(4)>=87 & x(4)<=87.5,interval(30,35)) ; ...
    };
benches{end+1} = bench;

bench.kfModel = @modelF16;
x = stl('x',16);
bench.requirements = {; ...
    %     "F16", globally(x(12)>0,interval(0,15)) ; ...
    };
benches{end+1} = bench;

bench.kfModel = @modelPowertrain;
x = stl('x',1);
u = stl('u',2);
rise = (u(1) < 8.8) & finally(u(1) > 40.0,interval(0,0.05));
fall = (u(1) > 40.0) & finally(u(1) < 8.8,interval(0,0.05));

rise2 = (u(1) < 8.8) & ~globally((u(1) < 40.0),interval(0,0.05));
fall2 = (u(1) > 40.0) & ~globally(u(1) > 8.8,interval(0,0.05));

beta=0.008;
gamma=0.007;
bench.requirements = {; ...
    "AFC_normal","(AFC27 0.008)",1, globally(implies(rise|fall,globally(abs(x(1))<beta,interval(1,5))),interval(11,50)); ...
    "AFC_normal","(AFC29 0.007)",1,globally(abs(x(1))<gamma,interval(11,50)) ; ...
    %     "AFC33",1, globally(abs(x(1))<gamma,interval(11,50)) ; ...
    };
benches{end+1} = bench;

for b = 1:length(benches)
    bench = benches{b};
    req = bench.requirements;
    for i = 1:size(req, 1)
        % initialize seeds
        rng(0)
        pyrunfile("seed.py")
        %         disp("--------------------------------------------------------")
        kfModel = bench.kfModel();
        property = req{i, 1};
        name = req{i, 2};
        kfModel.ak.dt = req{i, 3};
        eq = req{i, 4};

        if name == "NNx"
            kfModel.U = interval(1.95,2.05);
            kfModel.T=3;
        end
        if name == "AFC33"
            kfModel.U = interval([61.2;900],[81.2;1100]);
        end
        kfModel.spec = specification(eq,'logic');
        kfModel.runs=1;
        kfModel.verb=-1;
        kfModel.nResets=5; %'auto'
        kfModel.resetStrat=0;
        kfModel.maxSims=5000;

        for j=1:10
            [kfSolns,~] = falsify(kfModel);
            soln=kfSolns{1};
            fprintf('"%s","%s",%d,%d,%d,%g,%d,%g,%g',property,name,2,soln.sims,kfModel.T,soln.best.rob,soln.falsified,soln.simTime,soln.runtime);
            u=[soln.best.t,soln.best.u];
            [mrows, ncols] = size(u);
            fmturow=repmat('%g ',1,ncols);
            fmtucol=repmat(strcat(fmturow,';'),1,mrows);
            fprintf(strcat(',"[',fmtucol,']",'),u');
            fmtx0=repmat('%g ',1,numel(soln.best.x(1,:)));
            fprintf(strcat('"[',fmtx0,']" \n'),soln.best.x(1,:));
        end

        %print info
        %         fprintf('Benchmark: %s\n', name);
        %         printInfo(kfSolns)
    end
end
end