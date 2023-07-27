function Sys = KoopmanSetupAlpha(Sys)

%% setup
L=Sys.L; % horizon (# of steps)
cpBool=Sys.cpBool; %boolean of control points
reachZonos=Sys.reachZonos;

%% System dimensions and variables
nx=Sys.nx; %number of states
% variables
Sys.x = sdpvar(nx, L+1); %states
alpha = sdpvar(1, size(reachZonos{end}.generators,2));
if ~isempty(Sys.U)
    alphaU = alpha(size(reachZonos{1}.generators,2)+1:end);
    U = zonotope(Sys.U); c_u = center(U); G_u = generators(U);
    alphaU = reshape(alphaU,[size(G_u,2),length(alphaU)/size(G_u,2)]);
    c_u_ = repmat(c_u,1,size(alphaU,2));

    %append empty sdpvar for consistent length with states X
    addAlphaU = sdpvar(1,1);
    Falpha= -1<=addAlphaU<= 1;
    Sys.u = [c_u_ + G_u*alphaU, c_u+G_u*addAlphaU];
end

%constraints for alpha
Falpha= [Falpha,-1<=alpha<=1];

%constraint for control points
cpBool = cpBool(1:L,:); %get cpbool corresponding to number of steps
if ~isempty(cpBool) %piecewise constant signal and not pulse.
    %skip alphas that correspond to initial points
    k = size(reachZonos{1}.generators,2)+1;
    for row=1:size(cpBool,1)
        for col=1:size(cpBool,2)
            %if bool is zero constrain alpha to be same value as prev
            if ~cpBool(row,col)
                Falpha=[Falpha, alpha(k)==alpha(k-size(cpBool,2))];
            end
            k=k+1; %next alpha
        end
    end
end

%assign optim variables and outputs to system
Sys.Finit=Falpha; Sys.alpha=alpha;




