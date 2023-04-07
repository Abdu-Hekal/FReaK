function [tout, yout] = run_simulink(model_name, T, dt, x0, u)
    open_system(model_name,'loadonly');
    sim_mode = get_param(model_name,'SimulationMode');
    set_param(model_name,'SimulationMode','normal'); %need this to extract init state info, otherwise error is thrown if sim mode is not 'normal'

    xInitial = Simulink.BlockDiagram.getInitialState(model_name);
    for i=1:size(xInitial.signals,2)
        xInitial.signals(i).values = x0(i);
    end
    paramStruct.LoadInitialState = 'on';
    paramStruct.InitialState = 'xInitial';
    assignin('base','T',T);
    assignin('base','dt',dt);
    assignin('base','xInitial',xInitial);

    if ~isempty(find_system(model_name,'BlockType','Inport')) %check if simulink model has inputs
        paramStruct.LoadExternalInput= 'on';
        paramStruct.ExternalInput= 'u';
        assignin('base','u',u);
    end

    paramStruct.StopTime = 'T';
    paramStruct.FixedStep = 'dt';
    paramStruct.SaveTime='on';
    paramStruct.TimeSaveName='tout';
    paramStruct.SaveOutput='on';
    paramStruct.OutputSaveName='yout';
    paramStruct.SaveFormat='Array';
    paramStruct.SimulationMode=sim_mode;

    result = sim(model_name, paramStruct);
    tout = result.tout;
    yout = result.yout;

end
