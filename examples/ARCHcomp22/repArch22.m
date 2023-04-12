function solns = repArch22()
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
%     "AT2", globally(x(2) < 4750,interval(0,10)); ...
    %     "testAT2", globally(x(2) <= 4750,interval(0,10)); ...
        "AT51", globally(implies((x(3)>1) & finally((x(3)>=1 & x(3)<=1),interval(0.001,0.1)),finally(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),interval(0.001,0.1))),interval(0,30)); ...
    %      "AT6a", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
    %         "test", globally(x(1)<50 | x(1)>60,interval(10,30)),...
    %      "testAT6a", implies(globally(x(2)<3000,interval(0,4)),globally(x(1)<35,interval(0,4))); ...

    };

solns=dictionary(string.empty,cell.empty);
for i = 1:size(requirements, 1)
    for j = 1:2
        disp("--------------------------------------------------------")
        name = requirements{i, 1};
        eq = requirements{i, 2};

        model.spec = specification(eq,'logic');
        [model,~] = falsify(model);

        if model.soln.falsified
            disp(" ")
            fprintf("falsifying trace found! for requirement '%s'\n", name)
        else
            fprintf("No falsifying trace found! for requirement '%s'\n", name)
        end
        disp(['training iterations required: ',num2str(model.soln.trainIter)])

        if j==1
            solns(name)={{model.soln}};
        else
            soln=solns(name);
            soln{1}{end+1}=model.soln
            solns(name)=soln;
        end
    end
end
end

function avgTime = getAvgTime(soln)
    avgTime=0;
    for i=1:length(soln{1})
        avgTime=avgTime+soln{1}{i}.runtime;
    end
    avgTime=avgTime/i;
end

