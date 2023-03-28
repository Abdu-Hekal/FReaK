
x = stl('x',2);
eq = until(globally(x(2) <= 0.5,interval(0,0.2)),x(1)>=0.5,interval(0.4,1));
% eq = release(globally(x(2) < -0.2,interval(0,1)),x(1)>-0.9,interval(0,1));
% eq = finally(globally(x(2)>1.1,interval(0.4,0.6)),interval(0.4,1.2));

% eq = globally(finally(x(2)>-0.1,interval(0.1,0.3)),interval(0,1));
% eq = finally(x(2)>0.1,interval(0,1.2));

% eq = until(x(2) < 1,x(1)>=1.1,interval(0,1));


phi = -pi/2:0.01:0;
x = [cos(phi'),sin(phi')];
t = linspace(0,1,length(phi))';

% 
% phi = 0:0.1:1;
% x = [phi',phi'];
% t = linspace(0,1,length(phi))';


% eq = globally(implies((x(1)>0.5 & finally(x(2)>0.5,interval(0,0.1))),globally(x(2)>0.6, interval(0,0.2))),interval(0,1))



rho = computeRobustness(eq,x,t)
% checkStl(eq,x,t)
