function [tout, yout] = runSteamcondenser(T,~,u)
    
    assignin('base','u',u);
    assignin('base','T',T);
    
    result = sim('steamcondense_RNN_22', ...
        'LoadExternalInput', 'on', 'ExternalInput', 'u', ...
        'SaveTime', 'on', 'TimeSaveName', 'tout', ...
        'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
        'SaveFormat', 'Array');
    tout = result.tout;
    yout = result.yout;

end
