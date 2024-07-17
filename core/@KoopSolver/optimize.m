function Sys = optimize(Sys,options)
% OPTIMIZE Optimize the Koopman Solver (KoopSolver) object.
%
% Syntax:
%    Sys = optimize(Sys, options)
%
% Description:
%    This function optimizes the Koopman solver (KoopSolver) object to find 
%    falsifying trajectories. The optimization is done either directly or 
%    using an optimizer object. If no optimizer object is available, it
%    optimizes the system directly; otherwise, it uses the optimizer object.
%
% Inputs:
%    Sys - Koopman solver (KoopSolver) object
%    options - Options for the optimization process
%
% Outputs:
%    Sys - Updated Koopman solver (KoopSolver) object after optimization
%
% Example:
%    Sys = optimize(Sys, options);
%
% See also: setupOptimizer
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

if isempty(Sys.optimizer) %no optimizer object, optimize directly
    constraints=[Sys.Finit, Sys.Fstl, Sys.Fdyn];
    objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
    %% call solverarch
    optimize(constraints,objective,options);
else
    param = zeros(1,length(Sys.Ostl));
    if numEntries(Sys.offsetMap) > 0 %we have an offset
        keys = Sys.offsetMap.keys;
        param(keys)=Sys.offsetMap.values;
    end
    if isempty(Sys.Wstl)
        param={param};
    else
        %weighted stl is used, so pass weights also as param for solver
        %time steps
        timeIdxs = floor(Sys.solverTimePoints/Sys.koopdt)+1;
        weights=Sys.weights(:,timeIdxs);
        %make sure weights aren't too big
        weights=min(weights,Sys.maxWeight);
        param={param,weights};
    end
    [sol_control, errorflag1,~,~,P] = Sys.optimizer{param}; %% call solver
    assign(Sys.x,double(sol_control{1}));
    assign(Sys.Pstl,double(sol_control{2}));
    if ~isempty(Sys.reachZonos)
        assign(Sys.alpha,double(sol_control{3}));
        if ~isempty(Sys.u)
            assign(Sys.u,double(sol_control{4}));
        end
    else
        assign(Sys.u,double(sol_control{3}));
    end
end
end