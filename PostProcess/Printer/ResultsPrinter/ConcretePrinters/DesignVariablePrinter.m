classdef DesignVariablePrinter < ResultsPrinter
    
    properties (Access = protected, Abstract)
        simulationStr
        fieldName
    end
    
    properties (Access = protected)
        hasGaussData = false;
    end
    
    methods (Access = public)
        
        function printResults(obj,iter,fileID)
            f = obj.fields;
            dS = obj.createScalarDataBase(iter,fileID,f,obj.fieldName,'OnNodes');
            ScalarPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.x;
        end
        
    end
    
end