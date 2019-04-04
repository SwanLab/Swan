classdef NumericalHomogenizerDataBase < handle
    
    properties (SetAccess = private, GetAccess = public)
       dataBase                
    end
        
    properties (Access = private)
       femFileName
    end
    
    methods (Access = public)
        
        function obj = NumericalHomogenizerDataBase(femFileName)
            obj.femFileName = femFileName;
            obj.createDataBase();
        end
        
    end
    
    
    methods (Access = private)
        
        function createDataBase(obj)
            d = obj.createNumericalHomogenizerDataBase();
            d.levelSetDataBase       = obj.createLevelSetDataBase();
            d.interpDataBase         = obj.createInterpDataBase();
            d.volumeShFuncDataBase   = obj.createShVolumeDataBase();
            obj.dataBase = d;
        end
        
        function d = createNumericalHomogenizerDataBase(obj)
            edt = 'ElementalDensityCreatorByLevelSetCreator';
            d.elementDensityCreatorType = edt;
            d.outFileName  = 'NumericalHomogenizer';
            d.testName     = obj.femFileName;
            d.print = false;
            d.iter  = 0;
            d.pdim = '2D';
        end
        
        function d = createShVolumeDataBase(obj)
            d = SettingsShapeFunctional();
            d.filterParams.optimizer = 'MMA';
            d.filename = obj.femFileName;
            d.ptype = 'MICRO';
        end
        
    end
    
    methods (Access = private, Static)
        
        function d = createLevelSetDataBase()
            d = SettingsLevelSetCreator();
        end
        
        function d = createInterpDataBase()
            d = SettingsInterpolation();
        end
        
 
        
    end
    
    
end