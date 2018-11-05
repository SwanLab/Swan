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
        Density 
        direction
        LevelSetPostProcess
        DensityPostProcess
        Mesh
        OutPutName
        LevelOfNumFibers
        hasToBePrinted
        RotationMatrix
        DensityCreator
        densityPrinter
        angle
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
            obj.createDensity()
            obj.createDensityPrinter()
            obj.createMaterialProperties()
            obj.createRotationAngle()
            obj.setMaterialPropertiesInMicroProblem()
            obj.computeVolumeValue()
            obj.computeHomogenizedTensor()
            obj.print()
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
        
        function createDensity(obj)
            levFib = obj.LevelOfNumFibers;
            obj.DensityCreator = DensityCreatorByInitializer(levFib,obj.MicroProblem);
            obj.Density = obj.DensityCreator.getDensity();            
        end
        
        function createDensityPrinter(obj)
            quad = obj.MicroProblem.element.quadrature;
            mesh = obj.MicroProblem.mesh;
            obj.densityPrinter = DensityPrinter(quad,mesh);            
        end
        
        function createMaterialProperties(obj)
            obj.MatProps= obj.Interpolation.computeMatProp(obj.Density);
        end               
        
        function createRotationAngle(obj)
            fiberDirection = obj.direction.getValue;
            obj.angle = -acos(dot(fiberDirection,[1 0 0]));
        end
        
        function setMaterialPropertiesInMicroProblem(obj)
            obj.MicroProblem.setMatProps(obj.MatProps);
        end
        
        function computeVolumeValue(obj)
            volumeComputer = ShFunc_Volume(obj.Setting);
            vol = volumeComputer.computeCost(obj.Density);
            obj.Volume = vol;            
        end
        
        function computeHomogenizedTensor(obj)
           obj.MicroProblem.computeChomog();
           obj.rotateCh();
        end
        
        function rotateCh(obj)
            d = [ 0 0 1];
            dir = Vector3D;
            dir.setValue(d);
            a = obj.angle;
            Cv = obj.MicroProblem.variables.Chomog;
            C = SymmetricFourthOrderPlaneStressVoigtTensor();
            C.setValue(Cv);            
            obj.Ch = Rotator.rotate(C,a,dir);
        end
        
        function print(obj)
            if obj.hasToBePrinted
                d = obj.Density;
                outn = obj.OutPutName;
                obj.densityPrinter.print(d,outn)
            end
        end
    end
    
end

