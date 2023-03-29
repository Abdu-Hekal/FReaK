% Defining the plant dynamics
% The toolbox is organized around one main class, called STLClti. An
% STLClti object is primarily a continuous Linear Time Invariant (LTI)
% system with inputs, outputs and disturbances. Hence, a constructors for this
% class takes matrices to define such an LTI. We first define and A and B
% matrices for state evolution:

load("autokoopman_model.mat", "A","B")
Bu = B;
g = @(x) autokoopman(x);

R0=interval([0;1000;1],[0;1000;1]); 
U = interval([0;0],[100;325]); 

tic
R = reachKoopman(A,B,g,R0,U,1,0.01);
toc

%%
% Now we can call the main constructor of STLC_lti class.
Sys= Koopman_lti(R.zono);
%%
% In the next section, we will define the different settings for the
% control synthesis experiment. Before that, we define some initial state:

%% Defining the controller
% We start by defining the time instants for the whole experiment, the discrete time
% step ts for the controller and the horizon L in number of time steps.
Sys.time = 0:.01:1; % time for the dynamics
Sys.ts=0.01; % sampling time for controller
Sys.L=100;  % horizon (# of control inputs in blueSTL defined internally as (2*L)-1, AH modified so # of control inputs=L)

Sys.plot_x=2;
Sys.plot=true;


%%
% Then the following define a signal temporal logic formula to be satisfied
% by the system. Note that times in the temporal operators are continuous,
% not discrete steps.
% Sys.stl_list = {'ev_[0,1] x1(t) > 5', 'ev_[0,0.5] x2(t) > 5'};
% Sys.stl_list = {'ev_[0,1] x1(t) < 0.5 or ev_[0,1] x2(t) > 0.5'};
% Sys.stl_list = {'ev_[0,1] x1(t) < 0.5','ev_[0,1] x2(t) > 0.5'};
% Sys.stl_list = {'alw_[0,20] x1(t) < 120'};
% Sys.stl_list = {'alw_[0,10] x2(t) < 4750'};
Sys.stl_list = {'ev_[0.5,1] x2(t) < 1000 or ev_[0.5,1] x2(t) > 2000'};


%modify solver options:
solver = 'gurobi';  % gurobi, cplex, glpk
timeLimit = 2000;
gapLimit = 0.01;
gapAbsLimit = 0.1;
solnLimit = Inf;
verb = 2;
Sys.solver_options = sdpsettings('verbose', verb,'solver', solver, ...
    'gurobi.TimeLimit', timeLimit, ...
    'gurobi.MIPGap', gapLimit, ...
    'gurobi.MIPGapAbs', gapAbsLimit, ...
    'gurobi.SolutionLimit', solnLimit,...
    'cachesolvers',1,...
    'gurobi.BarHomogeneous', 1,...
    'gurobi.ScaleFlag', 2, ...
    'gurobi.DualReductions', 0);
Sys.bigM=100000000;

% Now we are ready to compile the controller for our problem.
tic
controller = get_controller(Sys);
setup_time = toc;
tic
Sys.run_open_loop(controller);
solve_time = toc;
% Sys.model_data.rob
% Sys.model_data.alpha'
% Sys.model_data.X(1,1:end)

elapsed_time = toc;
fprintf('Setup time: %f seconds\n', setup_time);
fprintf('Solve time: %f seconds\n', solve_time);

% alphaU = Sys.model_data.alpha(2:end);
% 
% U = zonotope(U); c_u = center(U); G_u = generators(U);
% alphaU = reshape(alphaU,[size(G_u,2),length(alphaU)/size(G_u,2)]);
% 
% u = c_u + G_u*alphaU

function R = reachKoopman(A,B,g,R0,U,tFinal,dt)
% compute reachable set for Koopman linearized model

    % compute initial set using Taylor model arithmetic
    n = dim(R0);
    tay = taylm(R0);
    tay = g(tay);
    disp(tay)
    R0 = polyZonotope(tay);
    R0 = polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:));

    % compute reachable set
    t = 0:dt:ceil(tFinal/dt)*dt; U = zonotope(U);

    set = cell(length(t),1); time = set; zono = set;

    set{1} = R0; time{1} = interval(-dt/2,dt/2);

    for i = 1:length(t)-1
        % AH edit to check if system has external input
        if B
            set{i+1} = A*set{i} + B*zonotope(U);
        else
            set{i+1} = A*set{i};
        end
        time{i+1} = time{i} + dt;
    end

    % compute ouput set
    for i = 1:length(set)
        set{i} = project(set{i},1:n);
        zono{i} = zonotope(set{i});
    end

    R.set = set; R.time = time; R.zono = zono;
end

