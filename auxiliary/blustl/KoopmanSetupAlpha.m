function Sys = KoopmanSetupAlpha(Sys)

%% setup
L=Sys.L; % horizon (# of steps)
cpBool=Sys.cpBool; %boolean of control points
reachZonos=Sys.reachZonos;

%% System dimensions and variables
nx=Sys.nx; %number of states
% variables
X = sdpvar(nx, L+1); %states
Alpha = sdpvar(1, size(reachZonos{end}.generators,2));

%constraints for Alpha
Falpha= -1<=Alpha<=1;

%constraint for control points
cpBool = cpBool(1:L,:); %get cpbool corresponding to number of steps
if ~isempty(cpBool) %piecewise constant signal and not pulse.
    %skip alphas that correspond to initial points
    k = size(reachZonos{1}.generators,2)+1;
    for row=1:size(cpBool,1)
        for col=1:size(cpBool,2)
            %if bool is zero constrain alpha to be same value as prev
            if ~cpBool(row,col)
                Falpha=[Falpha, Alpha(k)==Alpha(k-size(cpBool,2))];
            end
            k=k+1; %next alpha
        end
    end
end

%assign optim variables and outputs to system
Sys.Falpha=Falpha; Sys.Alpha=Alpha; Sys.X=X;




