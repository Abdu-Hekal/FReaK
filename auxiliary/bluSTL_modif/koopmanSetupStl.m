function Sys = koopmanSetupStl(Sys)

%% STL formula
var = struct('X',Sys.X);

stlList= KoopmanParseStlLabels(Sys);
M = Sys.bigM;
normz = Sys.normz;
phi = STLformula('phi', stlList{1});

%     [Fphi, Pphi] = KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M);
%     [Fphi, Pphi] = orig_KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M);
[Fstl, Pstl, Ostl] = vector_KoopmanMilpRobust(phi, 1, Sys.L+1, Sys.dt, var,M, normz);

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;




