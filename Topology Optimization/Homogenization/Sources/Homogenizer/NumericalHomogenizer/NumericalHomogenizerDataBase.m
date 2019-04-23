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
            d.materialInterpDataBase = obj.createMaterialInterpDataBase();
            d.materialDataBase       = obj.createMaterialDataBase();
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
            d.filterParams.filterType = 'P1';
            s = SettingsDesignVariable();
            s.type = 'Density';
            d.filterParams.designVar = DesignVariable.create(s);% Density(s);
            d.filename = obj.femFileName;
            d.materialInterpolationParams.type = 'ByInterpolation';
            d.domainType = 'MICRO';
        end
        
    end
    
    methods (Access = private, Static)
        
        function d = createLevelSetDataBase()
            d.levelSetType = 'horizontalFibers';
            d.levFib = 3;
            d.volume = 0.5;
        end
        
        function d = createMaterialInterpDataBase()
            d.materialInterpolation = 'SIMPALL';
        end
        
        function d = createMaterialDataBase()
            d.materialType = 'ISOTROPIC';
            d.matProp.rho_plus = 1;
            d.matProp.rho_minus = 0;
            d.matProp.E_plus = 1;
            d.matProp.E_minus = 1e-3;
            d.matProp.nu_plus = 1/3;
            d.matProp.nu_minus = 1/3;
        end        
        
    end
    
    
end