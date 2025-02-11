classdef DesignVariablePrinter < ResultsPrinter
    
    properties (Access = protected, Abstract)
        fieldName
    end
    
    properties (Access = protected)
        hasGaussData = false;
    end
    
    methods (Access = public)
        
        function printResults(obj,iter,fileID)
            f = obj.fields;
            dS = obj.createScalarDataBase(iter,fileID,f,obj.fieldName,'OnNodes');
            ScalarNodalPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields;
        end
        
    end
    
end