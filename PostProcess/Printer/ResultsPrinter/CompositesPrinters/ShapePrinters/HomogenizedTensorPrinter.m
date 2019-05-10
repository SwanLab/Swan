classdef HomogenizedTensorPrinter < CompositeResultsPrinter
    
    properties (Access = private)
        nstre
    end
    
    methods (Access = public)
        
        function obj = HomogenizedTensorPrinter(d)
            obj.simulationStr = 'HomogenizedTensor';            
            obj.computeNstre(d.ndim);
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            microProblems = d.phyProblems{1};
            fields = microProblems.variables2print;%  variables.var2print;
            for istre = 1:obj.nstre
                di.fields = fields{istre};
                p = obj.printers{istre};
                p.storeFieldsToPrint(di);
                p.setStrVariablesMicroCase(istre)                                
            end
        end
        
        function createPrinters(obj,d)
            for istre = 1:obj.nstre
                p = ResultsPrinter.create('ElasticityMicro',d);
                obj.printers{istre} = p;
                
            end
        end
        
    end
    
    methods (Access = private)
        
        function computeNstre(obj,ndim)
            if ndim == 2
                obj.nstre = 3;
            elseif ndim == 3
                obj.nstre = 6;
            end
        end
    end
    
    
end