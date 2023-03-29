classdef Koopman_lti < STLC_lti

    % Koopman reachability properties
    properties
        reach_zonos = []
        plot = false %plot traces
    end

    methods
        % Constructor
        function Sys = Koopman_lti(varargin)
            switch nargin
                case 1
                    reach_zonos = varargin{1};
                    nx = size(reach_zonos{1}.center,1);
                    nu = size(reach_zonos{2}.generators,2) - size(reach_zonos{1}.generators,2);
                    varargin{1} = zeros(nx);
                    varargin{2} = zeros(nx,nu);
            end
            Sys = Sys@STLC_lti(varargin{:})

            switch nargin
                case 1
                    Sys.reach_zonos = reach_zonos;
                    assert(reach_zonos{1}.isInterval, "initial set is not an interval")
                    Sys.x0 = [reach_zonos{1}.interval.inf,reach_zonos{1}.interval.sup];

            end
        end

        function controller = get_controller(Sys,enc)
            if nargin < 2
                if isempty(Sys.encoding)
                    enc = 'robust';
                else
                    enc = Sys.encoding;
                end
            end
            Sys.sysd = c2d(Sys.sys, Sys.ts);
            controller = reach_get_controller(Sys,enc);
        end

        function [Sys, status] = compute_input(Sys, controller)
            % computes the next input and update model data
            % status is 0 if everything is OK
            [Sys, status] = Koopman_compute_input(Sys, controller);
        end

        % Executes the controller in open loop mode
        function [Sys] = run_open_loop(Sys, controller)

            Sys = Sys.reset_data();
            rfprintf_reset();
            [Sys, status] = Sys.compute_input(controller);

            if status==0 && Sys.plot
                for v = 1:length(Sys.plot_x)
                    plot(Sys.model_data.time,Sys.model_data.X(Sys.plot_x(v),:));
                end
            end
        end
        function Sys = apply_input(Sys)
            Sys = Koopman_apply_input(Sys);
        end

    end
end


