classdef NewBoundaryConditions < handle
    % The idea is to merge BoundaryConditions and BCApplier into a single
    % class for now.
    
    properties (Access = private)
    end
    
    methods (Access = public)

        function obj = NewBoundaryConditions(cParams)
            obj.init(cParams);
        end
        
    end

    methods (Access = private)

        function init(obj, cParams)
        end

    end

end

