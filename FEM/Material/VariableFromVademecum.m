classdef VariableFromVademecum < handle
    
    
    properties (Access = protected, Abstract)
        fieldName
    end
    
    properties (Access = protected)
        vadVariables
        values
        nPoints
        nParams
    end
    
    methods (Access = protected)
        
        function obj = VariableFromVademecum(cParams)
            obj.init(cParams);
           
        end
        

        

        

        
    end
    

end