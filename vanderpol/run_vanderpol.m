function [tout, yout] = run_vanderpol(init_set, T)
    open_system('vanderpol','loadonly');
    set_param('vanderpol','SaveFormat','Structure');
    xInitial = Simulink.BlockDiagram.getInitialState('vanderpol');
    xInitial.signals(1).values = init_set(1);
    xInitial.signals(2).values = init_set(2);
            
    assignin('base','xInitial',xInitial);
    assignin('base','T',T);
    
    result = sim('vanderpol', ...
        'StopTime', 'T', ...
        'LoadInitialState', 'on', 'InitialState', 'xInitial', ...
        'SaveTime', 'on', 'TimeSaveName', 'tout', ...
        'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
        'SaveFormat', 'Array');
    tout = result.tout;
    yout = result.yout;
end
