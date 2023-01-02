classdef ScalarGaussVariablePrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr
        hasGaussData = true;
    end
    
    properties (Access = private)
        fieldName
    end
    
    methods (Access = public)
        
        function obj = ScalarGaussVariablePrinter(d)
            obj.init(d);
            obj.fieldName = d.fieldName;
            obj.simulationStr = [d.fieldName,'Gauss'];
        end
        
        function printResults(obj,iter,fileID)
            field(:,1,:) = squeeze(obj.fields)';
            dS = obj.createScalarGaussDataBase(iter,fileID,field, obj.fieldName,'OnGaussPoints');
            ScalarGaussPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields;
        end
        
        function createHeadPrinter(obj,d,dh)
            obj.headPrinter = GaussHeadPrinter(d,dh);
        end
        
    end
    
end