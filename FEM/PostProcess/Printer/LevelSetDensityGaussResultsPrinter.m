classdef LevelSetDensityGaussResultsPrinter < ResultsPrinter    
    properties
        fieldNameDensity = 'RegularizedDensity';
        fieldNameLevelSet = 'LevelSet';
        simulationCase = 'LevelSetDensityGauss'
        headPrinter = GaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = LevelSetDensityGaussResultsPrinter()
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
            iS = obj.istep;
            gaussDescriptor = 'Guass up?';
            dD = obj.createScalarGaussDataBase(obj.fileID,dens, obj.fieldNameDensity,iS,'OnGaussPoints',gaussDescriptor);
            dL = obj.createScalarDataBase(obj.fileID,ls, obj.fieldNameLevelSet,iS,'OnNodes');
            ScalarPrinter(dL);
            ScalarGaussPrinter(dD);  
            
        end
    
    end
    
    methods (Access = private)
        
        function d = createScalarGaussDataBase(obj,fileID,fieldValues,fieldName,istep,fieldPosition,gaussDescriptor)
            d.fileID = fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.gaussDescriptor = gaussDescriptor;
            d.simulationCase = obj.simulationCase;
        end        
        
        function d = createScalarDataBase(obj,fileID,fieldValues,fieldName,istep,fieldPosition)
            d.fileID = fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.simulationCase = obj.simulationCase;
        end
        
        
    end
    
    
end