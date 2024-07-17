function disp(obj)
% disp - override disp function to show object on command window
%
% Syntax:  
%    disp(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    ---
%
% Example: 
%    x = stl('x',2)
%    until(x(1) < 5,x(2) < 3,interval(0.1,0.2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Abdelrahman Hekal
% Written:      13-November-2023 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    disp(str(obj));
end