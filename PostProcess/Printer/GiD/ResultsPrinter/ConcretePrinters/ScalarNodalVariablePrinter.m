classdef ScalarNodalVariablePrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr
        hasGaussData = false;
    end
    
    properties (Access = private)
        fieldName
    end
    
    methods (Access = public)
        
        function obj = ScalarNodalVariablePrinter(d)
            obj.init(d);
            obj.fieldName = d.fieldName;
            obj.simulationStr = d.fieldName;
        end
        
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