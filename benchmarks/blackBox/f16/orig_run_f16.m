function [tout, yout] = orig_run_f16(altg, Vtg, phig, thetag, psig, T)
    model_err = false;
    analysisOn = false;
    printOn = false;
    plotOn = false;
    backCalculateBestSamp = false;

    powg = 9;                   % Power
    % Default alpha & beta
    alphag = deg2rad(2.1215);   % Trim Angle of Attack (rad)
    betag = 0;                  % Side slip angle (rad)

    t_vec = 0:0.01:T;

    % Set Flight & Ctrl Limits (for pass-fail conditions)
    [flightLimits,ctrlLimits,autopilot] = getDefaultSettings();
    ctrlLimits.ThrottleMax = 0.7;   % Limit to Mil power (no afterburner)
    autopilot.simpleGCAS = true;    % Run GCAS simulation

    % Build Initial Condition Vectors
    initialState = [Vtg alphag betag phig thetag psig 0 0 0 0 0 altg powg];
    orient = 4;             % Orientation for trim

    % Select Desired F-16 Plant
    % Table Lookup
    initialState
    [output, passFail] = RunF16Sim(initialState, t_vec, orient, 'stevens',...
        flightLimits, ctrlLimits, autopilot, printOn, plotOn);

	tout = t_vec;
	yout = output(12,:);
    vpa(output')
end
