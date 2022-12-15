classdef SectionVariablesComputer < handle
    
    properties (GetAccess = public, SetAccess = protected)
        designVariable
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public, Static)
        
        function section = create(cParams)
            f = SectionVariableFactory();
            section = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
          
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end
        
    end
    
end