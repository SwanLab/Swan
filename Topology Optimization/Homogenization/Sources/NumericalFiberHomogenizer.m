classdef NumericalFiberHomogenizer < handle
    
    properties
        Ch        
        MaterialValues
        Volume
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
        direction
        LevelSetPostProcess
        DensityPostProcess
        Mesh
        OutPutName
        LevelOfNumFibers
        hasToBePrinted
    end
    
    methods
        
        function obj = NumericalFiberHomogenizer(...
                      direction, LevelOfFibers,OutPutName, hasToBePrinted )
            obj.direction = direction;
            obj.OutPutName = OutPutName;
            obj.hasToBePrinted = hasToBePrinted;
            obj.LevelOfNumFibers = LevelOfFibers;
            obj.loadFileName();
            obj.loadTestData();
            obj.loadDimension();
            obj.createInterpolation();
            obj.createMicroProblem();
            obj.createLevelSet();
            obj.createFilter()
            obj.createDensity()
            obj.createMaterialProperties()
            obj.createPostProcess()
            obj.printLevelSet()
            obj.printDensity()
            obj.setMaterialPropertiesInMicroProblem()
            obj.computeVolumeValue()
            obj.computeHomogenizedTensor()
        end
    end
    
    methods (Access = private)
        
        function loadFileName(obj)
            obj.FileName = 'test_microFineFine';
        end
        
        function loadTestData(obj)
            obj.FileNameWithPath = strcat('./Input/',obj.FileName);            
            obj.Setting = Settings(obj.FileNameWithPath);
        end
        
        function loadDimension(obj)
            obj.Setting.pdim = '2D';
        end
        
        function createInterpolation(obj)
            MatValues       = obj.Setting.TOL;
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
            obj.Mesh = obj.MicroProblem.mesh;
            epsilon = obj.MicroProblem.mesh.mean_cell_size;
            LS_initializer = DesignVaribleInitializer_orientedFiber(...
                             obj.Setting,obj.Mesh,epsilon,...
                             obj.direction,obj.LevelOfNumFibers);            
            LS_initializer.compute_initial_design();
            obj.LevelSet = LS_initializer.x;
        end
        
        function createFilter(obj)
             dim = obj.Setting.ptype;
             fileName = obj.Setting.filename;
             obj.Filter = Filter_P1_LevelSet_2D(fileName,dim);
             obj.Filter.preProcess();            
        end
        
        function createDensity(obj)
            obj.Density = obj.Filter.getP0fromP1(obj.LevelSet);
        end
        
        function createMaterialProperties(obj)
            obj.MatProps= obj.Interpolation.computeMatProp(obj.Density);
        end
        
        function createPostProcess(obj)
            obj.createLevelSetPostProcess();
            obj.createDensityPostProcess();
        end
        
        function createLevelSetPostProcess(obj)
            SetOpt = 'SLERP';
            obj.LevelSetPostProcess = Postprocess_TopOpt.Create(SetOpt);
        end
        
        function createDensityPostProcess(obj)
            Quadrature = obj.MicroProblem.element.quadrature; 
            obj.DensityPostProcess = PostprocessDensityInGaussPoints(Quadrature);   
        end
        
        function printLevelSet(obj)
            if obj.hasToBePrinted
                results.iter = 0;
                results.case_file = obj.OutPutName;
                results.design_variable = obj.LevelSet;
                obj.LevelSetPostProcess.print(obj.Mesh,results);
            end
        end
        
        function printDensity(obj)
          if obj.hasToBePrinted
             results.iter = 0;             
             results.case_file = strcat(obj.OutPutName,'Density');
             results.density = obj.Density;
             obj.DensityPostProcess.print(obj.Mesh,results)          
          end
        end
        
        function setMaterialPropertiesInMicroProblem(obj)
            obj.MicroProblem.setMatProps(obj.MatProps);
        end
        
        function computeVolumeValue(obj)
            obj.Volume = obj.Filter.computeInteriorVolume(obj.LevelSet);
        end
        
        function computeHomogenizedTensor(obj)
           obj.MicroProblem.computeChomog();
           obj.Ch = obj.MicroProblem.variables.Chomog;
        end
        
    end
    
end

