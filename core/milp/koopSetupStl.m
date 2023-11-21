function Sys = koopSetupStl(Sys,hardcoded)


%% STL formula
n = Sys.solverdt/Sys.koopdt; %evaluate stl formula on x every n koopman time steps
x = Sys.x(:,1:n:end);
u = Sys.u(:,1:n:end);
var = struct('x',x,'u',u);
L=size(x,2);

phi= Sys.stl;
M = Sys.bigM;

if Sys.normalize
    normz = Sys.normz;
else
    normz = [];
end

if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
end
[Fstl, Pstl, Ostl] = koopStl(phi,1,L,Sys.solverdt,var,M,normz,hardcoded,Sys.offsetMap);

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;




