classdef NumericalFiberHomogenizer < NumericalHomogenizer
    
    
    properties (Access = private)      
        direction        
        levelOfNumFibers        
    end
    
    methods (Access = public)
        
        function obj = NumericalFiberHomogenizer(dir,level,fileName,print,iter)
            obj.init(fileName,print,iter);            
            obj.saveInputValues(dir,level);
            obj.generateMicroProblem();
            obj.computeHomogenizedVariables();
            obj.rotateCh();
            obj.createDensityPrinter();            
            obj.print()
        end        

    end
    
    methods (Access = protected)
        
        function createDensity(obj)
            levFib = obj.levelOfNumFibers;
            densityCreator = DensityCreatorByInitializer(levFib,obj.microProblem,obj.setting);
            obj.density = densityCreator.getDensity();
            obj.levelSet = densityCreator.getLevelSet();
        end       
        
    end
    
    methods (Access = private)
        
        function saveInputValues(obj,dir,levelFib)
            obj.direction = dir;
            obj.levelOfNumFibers = levelFib;
        end                                 
        
        function rotateCh(obj)
            angle = obj.computeRotationAngle();
            dir   = obj.computeZdirection();
            C     = obj.obtainChTensor();
            obj.Ch = Rotator.rotate(C,angle,dir);
        end
        
        function angle = computeRotationAngle(obj)
            fiberDirection = obj.direction.getValue;
            angle = -acos(dot(fiberDirection,[1 0 0]));
        end
        
        function C = obtainChTensor(obj)
            Cv = obj.Ch;
            C = SymmetricFourthOrderPlaneStressVoigtTensor();
            C.setValue(Cv);
        end
        
    end
    
    methods (Access = private, Static)
        
        function dir = computeZdirection()
            d = [ 0 0 1];
            dir = Vector3D;
            dir.setValue(d);     
        end        
      
    end
    
end

