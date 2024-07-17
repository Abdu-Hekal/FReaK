function [tout, yout] = runPowertrain(T,~,u)

    assignin('base','simTime',T);
    assignin('base','measureTime',1);
    assignin('base','fault_time',60);
    assignin('base','spec_num',1);
    assignin('base','fuel_inj_tol',1);
    assignin('base','MAF_sensor_tol',1);
    assignin('base','AF_sensor_tol',1);

    assignin('base','u',u);
    assignin('base','T',T);
    
    result = sim('AbstractFuelControl_M1', ...
        'StopTime', 'T', ...
        'LoadExternalInput', 'on', 'ExternalInput', 'u', ...
        'SaveTime', 'on', 'TimeSaveName', 'tout', ...
        'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
        'SaveFormat', 'Array');
    tout = result.tout;
    yout = result.yout(:,1);
end
