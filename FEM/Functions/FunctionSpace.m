classdef FunctionSpace < handle
    
    properties (Access = private)
        boundaryConditions
    end
    
    methods (Access = public)

        function obj = FunctionSpace(cParams)
            obj.init(cParams)
        end

        function assign(obj, bc)
            obj.boundaryConditions = bc;
        end

    end
    
    methods (Access = private)

        function init(obj, cParams)
            
        end
        
    end

end

