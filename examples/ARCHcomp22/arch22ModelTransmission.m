function results = arch22ModelTransmission()
% arch22ModelTransmission - runs all requirement formula for the  
%  model transmission benchmark of the ARCH'22 falsification Category
%
% Syntax:
%   results = arch22ModelTransmission()
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
max_train_size=100; %maximum number of training trajectories before quitting

x = stl('x',3);
requirements = {; ...
%     "AT1", globally(x(1) < 120,interval(0,20)),[1,2]; ...
    "AT2", globally(x(2) < 4750,interval(0,10)),[1,2]; ...
%     "AT51", globally(implies((x(3)>1) & next((x(3)>=1 & x(3)<=1),0),next(globally(x(3)>=1 & x(3)<=1,interval(0,2.5)),0)),interval(0,30)),[3,1]; ...
%      "AT6a", implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))),[1,2]; ...
%         "test", globally(x(1)<50 | x(1)>60,interval(10,30)),[1,2],...
%      "testAT6a", implies(globally(x(2)<3000,interval(0,4)),globally(x(1)<35,interval(0,4))),[1,2]; ...
    };

for i = 1:size(requirements, 1)
    disp("--------------------------------------------------------")
    
    name = requirements{i, 1};
    eq = requirements{i, 2};
    plot_vars = requirements{i, 3};

    model.spec = specification(eq,'logic');

    [falsified, trainset, crit_x, train_iter] = coreFalsify(model, max_train_size);

    if falsified
        disp(" ")
        fprintf("falsifying trace found! for requirement '%s'\n", name)
    else
        fprintf("No falsifying trace found! for requirement '%s'\n", name)
    end
    visualize_falsification(crit_x, trainset.t{1}, model.spec, plot_vars)
    disp(['training iterations: ',num2str(train_iter)])
    visualize_train(trainset, plot_vars)
    results = [falsified;train_iter];
end

