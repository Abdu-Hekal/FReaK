function Sys = koopmanSetupStl(Sys,hardcoded)

%% STL formula
var = struct('X',Sys.X);

stlList= KoopmanParseStlLabels(Sys);
M = Sys.bigM;
normz = Sys.normz;
phi = STLformula('phi', stlList{1});

%     [Fphi, Pphi] = KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M);
%     [Fphi, Pphi] = orig_KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M);
if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
    [Fstl, Pstl] = hardCodedvectorKoopmanMilpRobust(phi,1,Sys.L+1,Sys.dt,var,M,normz,Sys.offset,Sys.offsetCount);
    Ostl = {};
else
    [Fstl, Pstl, Ostl] = vector_KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M, normz);
end

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;




