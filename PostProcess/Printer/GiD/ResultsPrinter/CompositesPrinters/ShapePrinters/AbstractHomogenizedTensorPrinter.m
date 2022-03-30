classdef AbstractHomogenizedTensorPrinter < CompositeResultsPrinter
    
    properties (Access = protected)
        nstre
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.storeMicroProblemsFields(d);
          %  obj.storeRegularizedDensity(d);
        end
        
        function createPrinters(obj,d)
            obj.printers        = obj.createMicroProblemsPrinters(d);
          %  obj.printers{end+1} = obj.createRegularizedDensityPrinter(d);
        end
        
        function computeNstre(obj,ndim)
            if ndim == 2
                obj.nstre = 3;
            elseif ndim == 3
                obj.nstre = 6;
            end
        end
        
    end
    
    methods (Access = private)
        
        function p = createMicroProblemsPrinters(obj,d)
            for istre = 1:obj.nstre
                p{istre} = ResultsPrinter.create('ElasticityMicro',d);
                obj.printers{istre} = p;
            end
        end
        
        function p = createRegularizedDensityPrinter(obj,d)
            p = ResultsPrinter.create('DensityGauss',d);
        end
        
        function storeRegularizedDensity(obj,d)
            d.fields = d.regDensity;
            obj.printers{obj.nstre+1}.storeFieldsToPrint(d);
        end
        
    end
    
    methods (Access = protected, Abstract)
       storeMicroProblemsFields(obj)
    end
    
end