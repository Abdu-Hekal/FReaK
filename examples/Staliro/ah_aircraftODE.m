% This is a file by AH for using S-TALIRO with aircraft ODE
clear

disp(' ')
disp('aircraft example. ')
disp(' ')

model = @blackboxAircraft;

disp(' ')
disp(' The initial conditions defined as a hypercube:')
init_cond = [200 260;-10 10;120 150] %#ok<*NOPTS>

disp(' ')
disp(' The constraints on the input signals defined as a hypercube:')
input_range = [34386 53973; 0 16]
disp(' ')
disp(' The number of control points for each input signal:')
cp_array = [10 10];


preds(1).str = 'c';
preds(1).A = [1 0 0; -1 0 0];
preds(1).b = [260; -250];
preds(2).str = 'd';
preds(2).A = [1 0 0; -1 0 0];
preds(2).b = [240; -230];

preds(3).str = 'e';
preds(3).A = [0 0 -1];
preds(3).b = [0];

%define two stl spec outlined in paper
phis = {'(([]_[1,1.5] c) -> ([]_[3,4] !d))','[]_[0,4] e'};
%run with two different optimizers
optms = {'SA_Taliro','UR_Taliro'}

for ii=1:length(phis)
    for jj=1:length(optms)
        % initialize seed
        rng(0)
        %select stl spec
        phi=phis{ii};

        disp(' ')
        disp('Total Simulation time:')
        time = 4

        disp(' ')
        disp('Create an staliro_options object with the default options.')
        opt = staliro_options();

        disp(' ')
        disp('Change options:')
        opt.runs = 10;
        opt.spec_space = 'X';

        opt.optimization_solver = optms{jj};
        opt.optim_params.n_tests = 5000;

        opt.interpolationtype={'pconst'}
        opt.dispinfo=0;

        disp(' ')
        disp('Running S-TaLiRo ...')
        tic
        results = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);
        toc

        disp(' ')

        for i=1:numel(results.run)
            display(['Number of samples in Run = ',num2str(results.run(i).nTests)])
        end
    end
end
