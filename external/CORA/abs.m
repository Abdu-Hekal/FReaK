function res = abs(obj)
% abs - overloads the abs operator for stl objects.
% add this file to CORA stl class.
%
% Syntax:  
%    res = abs(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = abs(x(1))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    % check input arguments
    if ~isa(obj,'stl') ||  obj.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    % construct resulting stl object
    res = obj;
    
    res.type = 'abs';
    res.lhs = obj;
    res.rhs = [];
    res.from = [];
    res.to = [];
end

%------------- END OF CODE --------------