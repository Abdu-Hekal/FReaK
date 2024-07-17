function [tout, yout] = runModulator(T,x0,u)

open_system('modulator_3rd_order','loadonly');
sim_mode = get_param('modulator_3rd_order','SimulationMode');
set_param('modulator_3rd_order','SimulationMode','normal'); %need this to extract init state info, otherwise error is thrown if sim mode is not 'normal'

xInitial = Simulink.BlockDiagram.getInitialState('modulator_3rd_order');

xInitial.signals(1).values = x0(3);
xInitial.signals(2).values = x0(1);
xInitial.signals(3).values = x0(2);

paramStruct.LoadInitialState = 'on';
paramStruct.InitialState = 'xInitial';
assignin('base','T',T);
assignin('base','xInitial',xInitial);

if ~isempty(find_system('modulator_3rd_order','BlockType','Inport')) %check if simulink model has inputs
    paramStruct.LoadExternalInput= 'on';
    paramStruct.ExternalInput= 'u';
    assignin('base','u',u);
end

paramStruct.StopTime = 'T';
paramStruct.SaveTime='on';
paramStruct.TimeSaveName='tout';
paramStruct.SaveOutput='on';
paramStruct.OutputSaveName='yout';
paramStruct.SaveFormat='Array';
paramStruct.SimulationMode=sim_mode;

result = sim('modulator_3rd_order', paramStruct);
tout=result.tout;
yout = [result.x1,result.x2,result.x3];

end
