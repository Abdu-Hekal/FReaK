%set seed
rng(0)
pyrunfile("seed.py")

kf = modelF16();
[t, x, u, simTime] = sampleSimulation(kf);
tak = (0:kf.ak.dt:kf.T)'; %define autokoopman time points
xak = interp1(t,x,tak,kf.trajInterpolation); %define autokoopman trajectory points

[kf,trainset]=initialize(kf);
trainset.t{end+1} = tak;
trainset.X{end+1} = xak';
trainset.XU{end+1} = u(:,2:end)';

koopModel=learnKoopModel(kf,trainset);
%setup
A = koopModel.A;
B = koopModel.B;
g = koopModel.g;
R0=kf.R0;

dims=[26,21];
n = R0.dim; dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0,6,names);
tay = g(tay);
R0 = polyZonotope(tay);
R0=zonotope(R0);
hold on;
plot(R0,dims)

R0_=split(kf.R0,6);
R0 = R0_{1};
n = R0.dim; dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
R0=zonotope(R0);
plot(R0,dims,'g')

R0_=split(kf.R0,6);
R0 = R0_{2};
n = R0.dim; dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
R0=zonotope(R0);
plot(R0,dims,'r')




