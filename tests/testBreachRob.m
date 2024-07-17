Bdata = BreachTraceSystem({'temperature', 'humidity'});

func = -pi/2:0.01:0;
temperature = cos(func'); % in Celsius deg
humidity = sin(func'); % in percents
time = linspace(0,1,length(func))';

trace = [time temperature humidity]; % combine into a trace, column oriented
Bdata.AddTrace(trace);

% phi = STL_Formula('phi', '(alw_[0,0.2] humidity[t]<=0.5) until_[0.5,1] temperature[t]>=0.5');
% phi = STL_Formula('phi', '(alw_[0,0.2] humidity[t]<=0.5) until_[0.5,1] humidity[t]<=0.5');

% phi = STL_Formula('phi', 'ev (humidity[t]>0.8)');
% phi = STL_Formula('phi', 'alw_[0.4,0.6] (humidity[t]>1.1)');
% phi = STL_Formula('phi', 'alw((alw humidity[t]>0.8) until (temperature[t]>-0.9))');

% phi = STL_Formula('phi', '(alw humidity[t]<-0.4) until temperature[t]>0.9');

% This returns NaN 
% phi = STL_Formula('phi', 'ev_[0.2,1] (alw_[0.4,0.6] (humidity[t]>1.1))');


% phi = STL_Formula('phi', 'alw_[0,1] (alw_[0.1,0.3] humidity[t]>-0.1) ')

phi = STL_Formula('phi', 'alw_[0,0.2] humidity[t]<=0.5 and (alw_[0,0.2] humidity[t]<=0.5 or alw_[0,0.2] temperature[t]>=0.5)');

Rphi = BreachRequirement(phi);


robustness = @(Bdata) Rphi.Eval(Bdata);
rob=robustness(Bdata)




