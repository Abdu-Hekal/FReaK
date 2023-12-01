%set seed
rng(0)
pyrunfile("seed.py")

kf = modelF16();
[t, x, x0, u, kf] = randSimulation(kf);
tak = (0:kf.ak.dt:kf.T)'; %define autokoopman time points
xak = interp1(t,x,tak,kf.trajInterpolation); %define autokoopman trajectory points

[kf,trainset]=initialize(kf);
trainset.t{end+1} = tak;
trainset.X{end+1} = xak';
trainset.XU{end+1} = u(:,2:end)';

[kf, koopModel]=learnKoopModel(kf,trainset);
%setup
A = koopModel.A;
B = koopModel.B;
g = koopModel.g;
R0=kf.R0;

n = dim(R0); dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0,6,names);
tay = g(tay);
R0 = polyZonotope(tay);
hold on;
plot(R0,[21,22])

R0_=split(kf.R0,6);
R0 = R0_{1};
n = dim(R0); dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
plot(R0,[21,22],'g')

R0_=split(kf.R0,6);
R0 = R0_{2};
n = dim(R0); dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
plot(R0,[21,22],'r')


% w1 = width(interval(zonotope(zeros(size(R0.c)),R0.Grest)))
% w2 = width(interval(polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:))))



