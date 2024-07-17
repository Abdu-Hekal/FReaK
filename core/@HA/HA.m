classdef HA
    % HA - class representing a hybrid automaton object for setting up a hybrid
    % automaton to falsify
    %
    % Syntax:
    %    obj = HA(init,locs,guards)
    %
    % Inputs:
    %    init - initial location for automata (integar).
    %    locs - cell array of loc (location) structs, where each loc is of
    %    the following structure:
    %       loc.num - location number
    %       loc.dyn - dynamics for the location, can be provided as
    %       simulink model, blackbox function, OdeFcn, or ode
    %       loc.invs - invariants for the location, provided as a combination
    %       of linear equations, defined using CORA stl variables
    %
    %    jumps - cell of jumps/transitions, where each jump is a struct of
    %    the following structure
    %       jump.from - Number of "from" location
    %       jump.to - Number of "to" location
    %       jump.guards - guards for the jump, provided as a combination
    %       of linear equations, defined using CORA stl variables (optional)
    %       jump.reset - reset of variables associated with the jump,
    %       provided as a dictionary object with variables defined with 
    %       CORA stl object as keys and values are linear reset. reset can 
    %       be a combination of constants and variables
    %
    % Outputs:
    %    obj - generated hybrid automaton
    %
    %
    % See also: KF

    % Author:      Abdelrahman Hekal
    % Written:      13-January-2024
    % Last update:  ---
    % Last revision:---

    %------------- BEGIN CODE --------------
    properties
        init
        locs
        jumps
    end

    methods
        % Constructor
        function obj = HA(init,locs,jumps)
            % Validate input arguments
            assert(isnumeric(init) && isscalar(init) && isinteger(init), 'init must be a numeric, integer, scalar');
            assert(iscell(locs) && all(cellfun(@isstruct, locs)), 'locs must be a cell array of location structs');
            assert(iscell(jumps) && all(cellfun(@isstruct, jumps)), 'jumps must be a cell array of jump structs');

            % Validate the structure of location dynamics and invariants
            for i = 1:numel(locs)
                loc=locs{i};
                assert(isfield(loc, 'num'), 'Each location struct must have a field "num"');
                assert(isfield(loc, 'dyn'), 'Each location struct must have a field "dyn"');
                assert(isfield(loc, 'invs'), 'Each location struct must have a field "invs"');
            end

            % Validate the structure of jumps
            for i = 1:numel(jumps)
                jump = jumps{i};
                assert(isfield(jump, 'from') && isnumeric(jump.from) && isscalar(jump.from) && isinteger(init), 'Each jump struct must have a numeric scalar "from" location');
                assert(isfield(jump, 'to') && isnumeric(jump.to) && isscalar(jump.to) && isinteger(init), 'Each jump struct must have a numeric scalar "to" location');
                % Validate the structure of resets
                if isfield(jump, 'reset')
                    assert(isa(jump.reset,'dictionary'), 'Jump reset must be a dictionary object');
                    assert(isa(jump.reset.keys,'stl'), 'Jump reset keys must be state variables of type CORA stl')
                end
            end
            %assign fields
            obj.init=init;
            obj.locs=locs;
            obj.jumps=jumps;
        end
    end
end