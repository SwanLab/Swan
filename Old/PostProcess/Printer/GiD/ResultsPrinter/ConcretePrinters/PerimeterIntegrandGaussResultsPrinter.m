classdef PerimeterIntegrandGaussResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'PerimeterIntegrand';
        hasGaussData = true;
    end
    
    properties (Access = private)
        fieldName = 'PerimeterIntegrand';
    end
    
    methods (Access = public)
        
        function obj = PerimeterIntegrandGaussResultsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            dens = obj.fields;
            dS = obj.createScalarGaussDataBase(iter,fileID,dens, obj.fieldName,'OnGaussPoints');
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