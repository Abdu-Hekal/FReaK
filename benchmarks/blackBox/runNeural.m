function [tout, yout] = runNeural(T,~,u)

    u_ts=0.001;
    
    assignin('base','u',u);
    assignin('base','T',T);
    assignin('base','u_ts',u_ts);
    
    result = sim('narmamaglev_v1', ...
        'StopTime', 'T', ...
        'LoadExternalInput', 'on', 'ExternalInput', 'u', ...
        'SaveTime', 'on', 'TimeSaveName', 'tout', ...
        'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
        'SaveFormat', 'Array');
    tout = result.tout;
    yout = result.yout;
end
