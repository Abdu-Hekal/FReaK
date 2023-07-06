function Sys = koopmanSetupStl(Sys,hardcoded)

%% STL formula
n = Sys.solverdt/Sys.koopdt; %evaluate stl formula on x every n koopman time steps 
x = Sys.x(:,1:n:end);
u = Sys.u(:,1:n:end);
var = struct('x',x,'u',u);
L=size(x,2);

stl= KoopmanParseStlLabels(Sys);
M = Sys.bigM;
normz = Sys.normz;
phi = STLformula('phi', stl);

if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
    [Fstl, Pstl] = hardCodedvectorKoopmanMilpRobust(phi,1,L,Sys.solverdt,var,M,normz,Sys.offsetMap);
    Ostl = {};
else
    [Fstl, Pstl, Ostl] = vector_KoopmanMilpRobust(phi, 1, L, Sys.solverdt, var,M, normz);
end

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;




