function Sys = koopSetupDynamics(Sys)

L=Sys.L;

A=Sys.A; B=Sys.B;
x=Sys.x; u=Sys.u; 


%% Dynamics constraints
assert(all(Sys.X0.inf==Sys.X0.sup),'Initial region "R0" must be a point if reachability is not used')
F = x(:,1) == Sys.g(center(Sys.X0));
for k=2:L+1
    % x = Ax + Bu 
    F = [F, x(:,k) == A*x(:,k-1) + B*u(:,k-1)];
end

%% Dynamics constraints
Sys.Fdyn=F;

