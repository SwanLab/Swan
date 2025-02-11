classdef testHorizontalTensorRotatedVsSequentialLaminate < handle
    
    properties (Access = protected)
        C0
        C1
        lamDir
        lamPar
        theta
        lamTensor
        rotHorTensor
    end
    
    properties (Access = private)
        angle
        horLamDir
        rotDir 
        horTensor
    end
    
    methods (Access = protected)
        
        function computeTest(obj)
            obj.init()
            obj.computeHorizontalLaminate()
            obj.rotateHorizontalLaminate()
            obj.computeLaminateDirectly()
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createLaminateParameters()
            obj.createRotationParams()
            obj.createHorizontalDirection()
            obj.createLaminateDirection()
            obj.createMaterialTensors()
        end
        
        function createLaminateParameters(obj)
            obj.lamPar = 1;
            obj.theta = rand(1);
        end
        
        function createRotationParams(obj)
            obj.rotDir = Vector3D();
            obj.rotDir.setValue([0; 0 ;1]);
            obj.angle = pi/2*rand(1);
        end
        
        function createHorizontalDirection(obj)
            obj.horLamDir = Vector3D();
            obj.horLamDir.setValue([0 ;1 ;0]);
        end
        
        function createLaminateDirection(obj)
            d = obj.rotDir;
            a  = obj.angle;
            dir2Rotate = Vector3D();
            dir2Rotate.setValue(obj.horLamDir.getValue());
            rotatedDir = Rotator.rotate(dir2Rotate,a,d);
            obj.lamDir = rotatedDir;
        end
        
        function createMaterialTensors(obj)
            E1 = 1;
            E0 = 1e-3*E1;
            nu1 = 1/3;
            nu0 = 1/3;
            obj.C1 = IsotropicConstitutiveTensor(E1,nu1);
            obj.C0 = IsotropicConstitutiveTensor(E0,nu0);
        end
        
        function computeHorizontalLaminate(obj)
            c0        = obj.C0;
            c1        = obj.C1;
            dir{1}    = obj.horLamDir;
            m1        = obj.lamPar;
            frac      = obj.theta;
            lam = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,m1,frac);
            obj.horTensor = lam.getPlaneStressHomogenizedTensor();
        end
        
        function rotateHorizontalLaminate(obj)
            dir = obj.rotDir;
            alpha = (pi - obj.angle);
            Chor = obj.horTensor;
            obj.rotHorTensor = Rotator.rotate(Chor,alpha,dir);
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeLaminateDirectly(obj)
    end
    
end

