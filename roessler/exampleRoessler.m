% load Koopman linearized model
load('paramRoessler.mat');
B = zeros(size(A,1),1);
g = @(x) observablesRoessler(x);
dt = 0.05;
A = expm(A*dt);

% define initial and input set as well as time horizon
R0 = interval([-0.05;-8.45;-0.05],[0.05;-8.35;0.05]);
U = interval(0,0.1);
tFinal = 6;

% define specification
spec = specification(halfspace([0 -1 0],-6.124),'unsafeSet');
spec = specification(interval([-10;6;-10],[-6;8;10]),'unsafeSet');

x = stl('x',3);
eq = until(x(1) < 9,x(2) > -3,interval(0,4));
spec = specification(eq,'logic');

% falsification
[x0,u] = falsifyFixedModel(A,B,g,dt,spec,R0,U,tFinal);

% simulate the resulting trajectory on the original nonlinear system
f = @(x,u) [-(x(2)+x(3));
            x(1) + 0.2*x(2);
            0.2 + x(3)*(x(1) - 5.7)];

sysNonlin = nonlinearSys(f);

paramsSim.x0 = x0;
paramsSim.tFinal = size(u,2)*dt;
paramsSim.u = u;

[t,x] = simulate(sysNonlin,paramsSim);

param.R0 = R0;
param.U = U;
param.tFinal = tFinal;

simRes = simulateRandom(sysNonlin,param);

% visualization
figure; hold on; box on;
xlim([-10,10]); ylim([-10,10]);
if ~strcmp(spec.type,'logic')
    plot(spec.set,[1,2],'r');
end
plot(simRes,[1,2],'k');
plot(x(:,1),x(:,2),'g','LineWidth',1);