classdef LevelSetDensityGaussResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        fieldNameDensity  = 'RegularizedDensity';
        fieldNameLevelSet = 'LevelSet';
        simulationCase    = 'LevelSetDensityGauss';
        headPrinter = GaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = LevelSetDensityGaussResultsPrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            d.fileID = obj.fileID;
            d.gaussDescriptor = obj.gaussDescriptor;
            d.etype = obj.etype;
            d.ngaus = obj.ngaus;
            d.ndim  = obj.ndim;
            d.posgp = obj.posgp;
            obj.headPrinter.print(d);
        end
        
        function printResults(obj)
            dens = obj.fields.density;
            ls   = obj.fields.levelSet;
            dD = obj.createScalarGaussDataBase(dens, obj.fieldNameDensity,'OnGaussPoints');
            dL = obj.createScalarDataBase(ls, obj.fieldNameLevelSet,'OnNodes');
            ScalarPrinter(dL);
            ScalarGaussPrinter(dD);
        end
        
    end
    
    
end