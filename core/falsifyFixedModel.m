function [x0,u, kfModel] = falsifyFixedModel(A,B,g,kfModel)
% determine the most critical initial point and input trajectory of a
% Koopman linearized model for the given trajectory
%
% Input arguments:
%
%   -A:        system matrix for the Koopman linearized model
%   -B:        input matrix for the Koopman linearized model
%   -g:        function handle to the observable function for the Koopman
%              linearized model
%   -model:    koopman falsification model
%
% Output arguments:
%
%   -x0:       initial state for the most critical trajectory
%   -u:        piecewise constant inputs for the most critical trajectory

% compute reachable set for Koopman linearized model
R = reachKoopman(A,B,g,kfModel);
% determine most critical reachable set and specification
kfModel = critAlpha(R,kfModel);

% extract most critical initial state and input signal
[x0,u] = falsifyingTrajectory(kfModel);

% %modification to test (delete me)
drawu=u(:,2:end)';
x = g(x0);
for i = 1:size(drawu,2)
    x = [x, A*x(:,end) + B*drawu(:,i)];
end
figure; hold on; box on;
for i=1:size(drawu,2)
    plot(R.zono{i})
end
plot(x(1,1:400),x(2,1:400),'r','LineWidth',2);
drawnow

end

% Auxiliary Functions -----------------------------------------------------

