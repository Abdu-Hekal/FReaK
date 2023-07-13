 function [tout, yout] = run_f16(T,x0,~)
    
    altg = 4040;...
    Vtg = 540;... 
    %uncertain initial states
    phig = x0(4,:);
    thetag = x0(5,:);
    psig = x0(6,:);

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
    [output, passFail] = RunF16Sim(initialState, t_vec, orient, 'stevens',...
        flightLimits, ctrlLimits, autopilot, printOn, plotOn);

	tout = t_vec;
	yout = output';
end
