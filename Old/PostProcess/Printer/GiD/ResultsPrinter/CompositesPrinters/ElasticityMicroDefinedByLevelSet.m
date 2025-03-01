classdef ElasticityMicroDefinedByLevelSet < CompositeResultsPrinter
      
    methods (Access = public)
        
        function obj = ElasticityMicroDefinedByLevelSet(d)
            obj.simulationStr = 'ElasticityMicroDefinedByLevelSet';
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            obj.printers{1} = ResultsPrinter.create('DensityGauss',d);
            obj.printers{2} = ResultsPrinter.create('LevelSet',d);
            obj.printers{3} = ResultsPrinter.create('ElasticityMicro',d);
        end
        
        function storeFieldsToPrint(obj,d)
            obj.storeDensityGaussFields(d);
            obj.storeLevelSetFields(d);
            obj.storeMicroFields(d);
        end
        
    end
    
    methods (Access = protected)
        
        function storeDensityGaussFields(obj,d)
            di.fields = d.dens;
            obj.printers{1}.storeFieldsToPrint(di);
        end
        
        function storeLevelSetFields(obj,d)
            di.x = d.levelSet;
            obj.printers{2}.storeFieldsToPrint(di);
        end
        
        function storeMicroFields(obj,d)
            d.fields = d.microVar;
            obj.printers{3}.storeFieldsToPrint(d);
            obj.printers{3}.setStrVariablesMicroCase([])
        end
        
    end
    
    
end