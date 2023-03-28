classdef Koopman_lti < STLC_lti

    % Koopman reachability properties
    properties
        reach_zonos = []
        plot = false %plot traces
    end

    methods
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
                current_time =0;
                while (current_time < Sys.model_data.time(end))
                    out = sprintf('time:%g', current_time );
                    rfprintf(out);
                    Sys = Sys.apply_input();
                    Sys = Sys.update_plot();
                    drawnow;
                    current_time= Sys.system_data.time(end);
                end
                fprintf('\n');
            end
        end
        function Sys = apply_input(Sys)
            Sys = Koopman_apply_input(Sys);
        end

    end
end


