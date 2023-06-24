function Sys = koopmanSetupStl(Sys,hardcoded)

%% STL formula
n = Sys.solverdt/Sys.koopdt; %evaluate stl formula on x every n koopman time steps 
x = Sys.x(:,1:n:end);
nu = n/((size(Sys.x,2)-1)/size(Sys.u,2)); %evaluate stl formula on u every nu koopman time steps. Note that size of nu might be different than n for pulse input.
assert(floor(nu)==nu,'check if pulse input that isdivisible by trajectory steps')
u = Sys.u(:,1:nu:end);
var = struct('x',x,'u',u);
L=size(x,2);

stl= KoopmanParseStlLabels(Sys);

M = Sys.bigM;
normz = Sys.normz;
phi = STLformula('phi', stl);

if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
    [Fstl, Pstl] = hardCodedvectorKoopmanMilpRobust(phi,1,L,Sys.solverdt,var,M,normz,Sys.offset,Sys.offsetCount);
    Ostl = {};
else
    [Fstl, Pstl, Ostl] = vector_KoopmanMilpRobust(phi, 1, L, Sys.solverdt, var,M, normz);
end

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;




