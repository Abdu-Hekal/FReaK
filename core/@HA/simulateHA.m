function [tout, yout] = simulateHA(obj,T, x0, u)
% simulateHA - Simulate the hybrid automaton.
%
% Syntax:
%    [tout, yout] = simulateHA(obj,T, x0, u)
%
% Description:
%    This function simulates the hybrid automaton.
%
% Inputs:
%    obj - hybrid automaton object
%    T - total simulation time
%    x0  - Initial state vector for simulation
%    u   - Input vector for simulation
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%
% Example:
%    [tout, yout, simTime] = simulateHA(obj,T, x0, u);
%
% See also: HA
%
% Author:      Abdelrahman Hekal
% Written:     13-January-2024
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

% Validate input arguments
assert(isnumeric(T) && isscalar(T) && T>0, 'Time horizon (T) must be defined as a positive numeric')

% Initialize simulation variables
tstart=0;
tout = tstart;
yout = x0';
curLoc = obj.locs(obj.init);

while tstart < T   
    invs=curLoc.inv;
    guards={};
    for j=1:length(obj.jumps)
        jump=obj.jumps{j};
        if jump.from==curLoc.num
            guards{end+1}=jump.guards;
        end
    end
end
end

function [value,isterminal,direction] = events(t,y,guards,invs)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
end
