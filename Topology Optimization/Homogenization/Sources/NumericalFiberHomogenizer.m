classdef NumericalFiberHomogenizer < handle
    
    properties
        Ch        
        MaterialValues
        VolumeValue
    end
    
    properties (Access = private)
        FileName
        FileNameWithPath
        Setting
        MicroProblem
        Interpolation
        MatProps
        Filter
        LevelSet
        Density 
        Volume
        direction
    end
    
    methods
        
        function obj = NumericalFiberHomogenizer(direction )
            obj.direction = direction;
            obj.loadFileName();
            obj.loadTestData();
            obj.loadDimension();
            obj.createInterpolation();
            obj.createMicroProblem();
            obj.createLevelSet();
            obj.createFilter()
            obj.createDensity()
            obj.createMaterialProperties()
            obj.setMaterialPropertiesInMicroProblem()
            obj.createVolumeFunctional()
            obj.computeVolumeValue()
            obj.computeHomogenizedTensor()
        end
    end
    
    methods (Access = private)
        
        function loadFileName(obj)
            obj.FileName = 'test_microHorizontalFine';
        end
        
        function loadTestData(obj)
            obj.FileNameWithPath = strcat('./Input/',obj.FileName);            
            obj.Setting = Settings(obj.FileNameWithPath);
           % obj.Setting.initial_case = 'orientedFiber';
           % obj.Setting.initial_case = 'horizontal';
        end
        
        function loadDimension(obj)
            obj.Setting.pdim = '2D';
        end
        
        function createInterpolation(obj)
            MatValues  = obj.Setting.TOL;
            Material        = obj.Setting.material;
            InterpFunction  = obj.Setting.method;
            Dim             = obj.Setting.pdim;
            Interp          = Material_Interpolation.create(MatValues,Material,InterpFunction,Dim);
            obj.Interpolation = Interp;
            obj.MaterialValues = MatValues;
        end
        
        function createMicroProblem(obj)                        
            obj.MicroProblem = Elastic_Problem_Micro(obj.Setting.filename);
            obj.MicroProblem.preProcess();
        end
        
        function createLevelSet(obj)
            Mesh = obj.MicroProblem.mesh;
            epsilon = obj.MicroProblem.mesh.mean_cell_size;
            LS_initializer = DesignVaribleInitializer_orientedFiber(obj.Setting,Mesh,epsilon,obj.direction);            
            LS_initializer.compute_initial_design();
            obj.LevelSet = LS_initializer.x;
        end
        
        function createFilter(obj)
             dim = obj.Setting.ptype;
             fileName = obj.Setting.filename;
             obj.Filter = Filter_P1_LevelSet_2D_Interior(fileName,dim);
             obj.Filter.preProcess();            
        end
        
        function createDensity(obj)
            obj.Density = obj.Filter.getP0fromP1(obj.LevelSet);
        end
        
        function createMaterialProperties(obj)
            obj.MatProps= obj.Interpolation.computeMatProp(obj.Density);
        end
        
        function setMaterialPropertiesInMicroProblem(obj)
            obj.MicroProblem.setMatProps(obj.MatProps);
        end
        
        function createVolumeFunctional(obj)
            obj.Volume = ShFunc_Volume(obj.Setting);
            obj.Volume.filter.preProcess();
        end
        
        function computeVolumeValue(obj)
            obj.Volume.computeCostAndGradient(obj.LevelSet)            
            obj.VolumeValue = obj.Volume.value;
        end
        
        function computeHomogenizedTensor(obj)
           obj.MicroProblem.computeChomog();
           obj.Ch = obj.MicroProblem.variables.Chomog;
        end
        
    end
    
end

