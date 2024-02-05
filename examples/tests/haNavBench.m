
%initial set
init.loc = 13;
init.cube = [0.2 0.8; 3.2 3.8; -0.4 0.4; -0.4 0.4];
%initial point
init.h0.l0=init.loc;
init.h0.t0=0;
%random initial continous state
init.h0.x0=rand(size(init.cube, 1), 1) .* (init.cube(:, 2) - init.cube(:, 1)) + init.cube(:, 1);
A = [4 2 3 4; 3 6 5 6; 1 2 3 6; 2 2 1 1];
%create nav bench hybrid automata
model = navbench_hautomaton(0,init,A);
model.init=init;

%create KF object for nav bench HA
kf=KF(model);
kf.R0 = interval([13;0.2;3.2;-0.4;-0.4],[13;0.8;3.8;0.4;0.4]);

kf.T=25;
kf.dt = 0.01;
kf.ak.dt=0.1; %2.5

%first x is loc, and other 4 is continous states
x = stl('x',5);
eq = globally(~(x(1)>=11 & x(1)<=11),interval(0,25));
kf.spec = specification(eq,'logic');
kf.verb=3;

[solns,allDatas]=falsify(kf);
soln=solns{1};
allData=allDatas{1};

visualizeTrain(allData,kf.ak.dt,1)



% [tout, yout]=sampleSimulation(kf);
