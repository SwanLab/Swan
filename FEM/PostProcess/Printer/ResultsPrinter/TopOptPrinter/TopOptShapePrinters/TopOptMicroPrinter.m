classdef TopOptMicroPrinter < TopOptShapePrinter
    
    properties (Access = private)
        nstre
    end
    
    methods (Access = public)
        
        function obj = TopOptMicroPrinter(d)
            obj.computeNstre(d.ndim);            
            obj.createPrinters(d);            
        end
        
        function itHas = hasGaussData(obj)
            itHas = true;
        end
        
        function storeResultsInfo(obj,shape)
            for istre = 1:obj.nstre
                phyPr = shape.getPhysicalProblem();
                d = obj.obtainVariablesAndQuad(phyPr,istre);
                obj.printers{istre}.storeResultsInfo(d);
            end
        end
       
        function print(obj,istep)
            for istre = 1:obj.nstre
                p = obj.printers{istre};
                p.setStrVariablesMicroCase(istre)                
                p.printOnlyResults(istep);
            end
        end
        
    end
    
    methods (Access = protected,Static)
    
        function d = obtainVariablesAndQuad(phyPr,istre)
            d.quad = phyPr.element.quadrature;
            d.variables = phyPr.variables.var2print{istre};                                    
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
        
        function createPrinters(obj,d)
            obj.printers = cell(obj.nstre,1);
            for istre = 1:obj.nstre
                sh = ResultsPrinter.create('ElasticityMicro',d);
                obj.printers{istre} = sh;
            end
        end
        
    end
    
end