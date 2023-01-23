t__ = [0; 10; 20; 30];
u__ = [50 0; 100 0; 0 325; 0 325];

p = [];
u = [t__, u__];

T = 30;

[tout, yout] = run_transmission(p, u, T);