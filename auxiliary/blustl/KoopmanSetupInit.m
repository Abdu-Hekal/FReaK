function Sys = KoopmanSetupInit(Sys)

L=Sys.L;
Sys.x = sdpvar(Sys.nx+Sys.nObs, L+1); %states
Sys.u = sdpvar(Sys.nu, L); %inputs

U=Sys.U; u=Sys.u; 

%% Input constraints
assert(isequal(class(U),'interval'), "U must be an interval if not using reachability (in the version)")

F= repmat(U.inf,1,L)<=u<=repmat(U.sup,1,L);

Sys.Finit=F;

