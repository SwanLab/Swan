classdef TopOptResultPrinterWithGauss < TopOptResultsPrinter
    
    
    methods (Access = public)
        
        function obj = TopOptResultPrinterWithGauss(d,dT)
            obj.init(d,dT);
        end
    end
    
    methods (Access = protected)
                
        function createHeadPrinter(obj,d)
                obj.headPrinter = GaussHeadPrinter;
                obj.headerData  = obj.createHeadPrinterDataBase(d);
        end
        
        function createDataBaseForResultsPrinter(obj,d)
            dI.dStandard = d.dStandard;
            dI.dStandard.ndim           = obj.ndim;
            dI.dStandard.simulationStr  = obj.simulationStr;
            dI.dGauss.ngaus = d.dGauss.ngaus;
            dI.dGauss.posgp = d.dGauss.posgp;
            dI.dGauss.gaussDescriptor = d.dGauss.gaussDescriptor;
            obj.dResultsPrinter = dI;
        end
    end
    
    methods (Access = private)
        
        function hD = createHeadPrinterDataBase(obj,d)
            hD.etype = obj.etype;
            hD.ndim  = obj.ndim;
            hD.gaussDescriptor = d.dGauss.gaussDescriptor;
            hD.ngaus = d.dGauss.ngaus;
            hD.posgp = d.dGauss.posgp;
        end
        
    end
    
    
    
end