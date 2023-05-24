classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        dt %time_step

        stlList %stl list to falsify
        cpBool %boolean array representing cp points

        %milp sdpvars and constraints
        Falpha %constraints on alpha
        Fstl %constraints for stl
        normFstl % normalized constraints for stl
        Freach %constraints states according to reachable sets
        X %states optim var
        Alpha
        Pstl %robustness of stl

        %optimizer object
        milp
    end

    properties (Dependent)
        bigM %bigM value for milp
        xlabel %default label name bluSTL
        nx %number of state variables
        L %number of control points
        normz %normalization values based on boundaries of reachable sets
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reachZonos,dt)
            Sys.reachZonos = reachZonos;
            Sys.dt=dt;

        end

        function Sys = setupAlpha(Sys)
            Sys = KoopmanSetupAlpha(Sys);
        end

        function Sys = setupStl(Sys)
            Sys = koopmanSetupStl(Sys);
        end

        function Sys = setupReach(Sys)
            Sys = KoopmanSetupReach(Sys);
        end

        function Sys = normStl(Sys)
            F = Sys.Fstl
            normz = Sys.normz;
            normFstl=[];
            for i = 1:length(F)
                predTag = tag(F(i));
                if contains(predTag, 'pred')
                    idxs=split(predTag,',');
                    normVal=0;
                    for j=2:length(idxs)
                        idx=str2double(idxs(j));
                        normVal = normVal + normz(idx);
                    end
                    normFstl=[normFstl,[sdpvar(F(i))/normVal >= 0]:predTag];
                else
                    normFstl=[normFstl,F(i)]; 
                end
            end
            Sys.normFstl = normFstl;
        end

        function diagnostics = optimize(Sys,options)

            constraints=[Sys.Falpha, Sys.Fstl, Sys.Freach];
            objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
            %% call solverarch
            diagnostics = optimize(constraints,objective,options);
        end
        %getters for dependent properties
        function bigM = get.bigM(Sys)
            %find suitable bigM based on zonotope boundaries
            bigM=0;
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                if norm(zono,inf) > bigM
                    order = ceil(log10(norm(zono)));
                    bigM = 10^order;
                end
            end
        end
        function nx = get.nx(Sys)
            nx=size(Sys.reachZonos{1}.center,1);
        end
        function L=get.L(Sys)
            L=size(Sys.reachZonos,1)-1;

        end
        function xlabel=get.xlabel(Sys)
            % default label names
            xlabel = cell(1,Sys.nx);
            for iX = 1:Sys.nx
                xlabel{iX} = ['x' num2str(iX)];
            end
        end
        function normz = get.normz(Sys)
            nx=Sys.nx;
            minBound=inf(nx,1);
            maxBound=-inf(nx,1);
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                for k=1:nx
                    supFun=zeros(1,nx);
                    supFun(k)=-1;
                    minBound(k) = min(minBound(k),-supportFunc(zono,supFun));
                    supFun(k)=1;
                    maxBound(k) = max(maxBound(k),supportFunc(zono,supFun));
                end
            end
            normz = maxBound-minBound;
        end
    end
end

