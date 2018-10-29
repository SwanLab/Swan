classdef testHorizontalTensorRotatedVsSequentialLaminate < test
    
    properties (Access = protected)
        C0
        C1
        lamDir  
        lamPar
        theta
        lamTensor
    end
    
    properties (Access = private)
        angle
        horLamDir
        rotDir 
        horTensor
        rotHorTensor
    end    
    
    properties (Access = protected,Abstract)
       tol 
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
            obj.createTolerance()
        end
        
        function createLaminateParameters(obj)
            obj.lamPar = 1;
            obj.theta = rand(1);
        end
                       
        function createRotationParams(obj)
            obj.rotDir = [0 0 1];
            obj.angle = pi/2*rand(1);
        end
        
        function createHorizontalDirection(obj)
            obj.horLamDir(:,1) = [0 1 0];
        end
        
        function createLaminateDirection(obj)
            d = obj.rotDir;
            a  = obj.angle;
            dir2Rotate = obj.horLamDir;
            rotatedDir = Rotator.rotate(dir2Rotate,a,d);            
            obj.lamDir = rotatedDir;
        end
        
        function createMaterialTensors(obj)
            E1 = 1;
            E0 = 1e-3*E1;
            nu1 = 1/3;
            nu0 = 1/3;
            obj.C1 = IsotropicConstitutiveTensor3D(E1,nu1);
            obj.C0 = IsotropicConstitutiveTensor3D(E0,nu0);
        end
        
        function computeHorizontalLaminate(obj)
            c0        = obj.C0;
            c1        = obj.C1;
            dir(1,:)  = obj.horLamDir;
            m1        = obj.lamPar;
            frac      = obj.theta;
            lam = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,m1,frac);
            obj.horTensor = lam.getPlaneStressHomogenizedTensor();
        end
        
        function rotateHorizontalLaminate(obj)
            dir = obj.rotDir;
            alpha = (pi - obj.angle);    
            Chor = obj.horTensor;
            C = FourthOrderVoigtTensor();
            C.setValue(Chor);
            obj.rotHorTensor = Rotator.rotate(C,alpha,dir);
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            rotHor = obj.rotHorTensor;
            lTens  = obj.lamTensor;
            tolerance = obj.tol;
            hasPassed = norm(rotHor - lTens) < tolerance;
        end
    end
    
    methods (Abstract, Access = protected)
        computeLaminateDirectly(obj)  
        createTolerance(obj)
    end
    
end

