classdef DensityResultsPrinter < ResultsPrinter
    
    
    properties (Access = protected)
        simulationStr = 'NodalDensity';
    end
    
    properties (Access = private)
        fieldName   = 'Density';
        headPrinter = NoGaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = DensityResultsPrinter(d)
            obj.init(d);
            obj.printHeader();
        end
        
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.headPrinter.print(obj.fileID);
        end
        
        function printResults(obj)
            dens = obj.fields;
            dS = obj.createScalarDataBase(dens,obj.fieldName,'OnNodes');
            ScalarPrinter(dS);
        end
        
    end
    
    
end