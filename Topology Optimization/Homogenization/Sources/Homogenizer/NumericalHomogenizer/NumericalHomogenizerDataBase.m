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
            d.volumeShFuncDataBase   = obj.createShVolumeDataBase(d);
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
        
        function d = createShVolumeDataBase(obj,dI)

            
            d = SettingsShapeFunctional();
            d.filterParams.filterType = 'P1';
            s = SettingsDesignVariable();
            s.type = 'Density';
            s.initialCase = 'full';
            

            
            fileName = s.mesh;
            dF = FemInputReader_GiD().read(fileName);
            cParams.coord  = dF.mesh.coord;
            cParams.connec = dF.mesh.connec;
            meshT = Mesh.create(cParams);
            
            s.mesh = meshT;

            scalarPr.epsilon = 1e-3;
            scalarPr.mesh = meshT;
            s.scalarProductSettings = scalarPr;

            % (19/12/2023): The future idea will be to destroy
            % LevelSerCreator and use GeometricalFunction
            sLs        = s.creatorSettings;
            sLs.ndim   = s.mesh.ndim;
            sLs.coord  = s.mesh.coord;
            sLs.type   = s.initialCase;
            lsCreator  = LevelSetCreator.create(sLs);
            phi        = lsCreator.getValue();
            switch s.type
                case 'Density'
                    value = 1 - heaviside(phi);
                case 'LevelSet'
                    value = phi;
            end
            ss.fValues = value;
            ss.mesh    = s.mesh;
            ss.order   = 'P1';
            s.fun      = LagrangianFunction(ss);

            designVar = DesignVariable.create(s);% Density(s);
            d.femSettings.fileName = obj.femFileName;
            d.femSettings.scale = 'MICRO';
            mesh = designVar.mesh;
            
            sHomog.type                   = 'ByInterpolation';
            sHomog.interpolation          = dI.materialInterpDataBase.materialInterpolation;
            sHomog.dim                    = dI.pdim;
            sHomog.typeOfMaterial         = dI.materialDataBase.materialType;
            sHomog.constitutiveProperties = dI.materialDataBase.matProp;
            sHomog.nelem                  = size(mesh.coord,1);
            sHomog = SettingsHomogenizedVarComputer.create(sHomog);
            d.homogVarComputer = HomogenizedVarComputer.create(sHomog);
        end
        
    end
    
    methods (Access = private, Static)
        
        function d = createLevelSetDataBase()
            d.type = 'horizontalFibers';
            d.levelOfFibers = 3;
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