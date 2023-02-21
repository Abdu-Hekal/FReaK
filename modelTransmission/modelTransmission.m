% generate training data
X = {}; U = {}; t = {}; N = 100;

for i = 1:10 
    u = [linspace(0,30,N)',100*rand(N,1),320*rand(N,1)];
    [t{end+1}, x] = run_transmission([], u, 30);
    U{end+1} = repelem(u(:,2:end)',1,30); X{end+1} = x';
end
dt = t{1}(2) - t{1}(1);

% identify Koopman model
numFeat = 100;                      % number observables
rank = 5;                           % rank for DMD
l = 0.001;                          % lengthscale for kernel

[cost,g,A,B,C,c] = costFunFourier([numFeat;rank;l;1],X,U,dt,false,false);                                    

% save observables in function
path = fileparts(which(mfilename()));
xSym = sym('x',[size(X{1},1),1]);
matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'functionTransmission'));

% visualize the predictions for the identfied Koopman model
sys = linearSysDT(linearSys(A,B,c,C),dt);

figure; hold on; box on;

for r = 1:length(X)
    x = simulateKoopman(sys,@(x) functionTransmission(x),X{r}(:,1),U{r});
    plot(X{r}(1,:),X{r}(2,:),'r');
    plot(x(1,:),x(2,:),'b');
end

% reachability analysis (reachable set for all possible inputs)
params.R0 = zonotope(functionTransmission(X{1}(:,1)));
params.U = zonotope(interval([0;0],[100;320]));
params.tFinal = 30;

options.zonotopeOrder = 10000;

R = reach(sys,params,options);

% determine most critical input from reachable set
ind = []; val = -inf;

for i = 1:length(R.timePoint.set)
    val_ = supportFunc(R.timePoint.set{i},[0 1 0]);
    if val_ > val
        val = val_; ind = i;
    end
end

[~,~,alpha] = supportFunc(R.timePoint.set{ind},[0 1 0]);
alpha = reshape(alpha,[2,round(length(alpha)/2)]);
u = center(params.U) + generators(params.U) * alpha;

% run most critical input on the real system
u_ = [linspace(0,30,size(u,2))', u'];

[~, x] = run_transmission([], u_, 30);
% 
plot([0 90],[4750, 4750],'--r')
plot(x(:,1),x(:,2),'g');