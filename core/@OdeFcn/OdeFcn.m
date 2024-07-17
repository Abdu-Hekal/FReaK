classdef OdeFcn
    % OdeFcn - Class for simulating ODEs using a function handle and solver options
    %
    % Syntax:
    %    obj = OdeFcn(handle, solver, options)
    %
    % Inputs:
    %    handle - Function handle for the ODEs (e.g., @(t, x) odefun(t, x, u))
    %    solver - Solver type (e.g., 'ode45', 'ode23', etc.)
    %    options - Options of type odeset for the chosen ODE solver (optional)
    %
    % Outputs:
    %    obj - Object of the OdeFcn class
    %
    % Example:
    % odeFun = @(t, x) -x; % Example ODE function
    % solver = 'ode45';
    % options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    % obj = OdeFcn(odeFun, solver, options);
    %
    % See also: ode45, ode23, odeset
    %
    % Author:      Abdelrahman Hekal
    % Written:     14-January-2024
    % Last update: ---
    % Last revision: ---

    properties
        handle % Function handle for the ODEs
        solver % Solver type (e.g., 'ode45', 'ode23', etc.)
        options % Options for the chosen ODE solver (optional)
    end

    methods
        % Constructor
        function obj = OdeFcn(handle, solver, options)
            % Validate input arguments
            narginchk(1, 3);

            % Assign properties
            assert(isa(handle, 'function_handle'),'handle property for the OdeFcn must be a function handle');
            assert(nargin(handle)>=2&&nargin(handle)<=3,['ODE function handle must accept 2 or 3 ' ...
                'input arguments, t, x and (optional) u']);
            obj.handle = handle;

           % Set solver type; use provided value or default to 'ode45'
            if nargin >= 2
                obj.solver = solver;
            else
                obj.solver='ode45';
            end

            % If options is provided, assign it; otherwise, use empty options
            if nargin == 3
                assert(isa(options,'struct'),'options for OdeFcn must be a struct of odeset')
                obj.options = odeset(options);
            else
                obj.options = odeset;
            end
        end
    end
end
