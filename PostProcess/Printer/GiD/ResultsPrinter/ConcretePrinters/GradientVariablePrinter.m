classdef GradientVariablePrinter < ResultsPrinter
    
    properties (Access = protected)
        fieldName = 'Gradient';
    end
    
    properties (Access = protected)
        hasGaussData = false;
        simulationStr
    end
    
    methods (Access = public)
        
        function obj = GradientVariablePrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            f = obj.fields;
            dS = obj.createScalarDataBase(iter,fileID,f,obj.fieldName,'OnNodes');
            ScalarNodalPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.gradient;
        end
        
    end
    
end