classdef testCommutingVoigtHomog < handle
    
    properties (Access = private)
        theta
        direction
        stiffTensor
        weakTensor
        vhpTensor
        hvpTensor
    end
    
    properties (Access = public)
        tol = 1e-12;
    end
    
    methods (Access = public)
        
        function obj = testCommutingVoigtHomog()
            obj.init()
            obj.computeWeakStiffTensor()
            obj.computeHomogVoigtPlaneStressTensor()
            obj.computeVoigtHomogPlaneStressTensor()
        end
       
         function error = computeError(obj)
            c1 = obj.vhpTensor.getValue();
            c2 = obj.hvpTensor.getValue();
            error = norm(c2-c1)/norm(c1);
         end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.theta = rand(1);
            angle = rand(1)*pi/2;
            dir = [cos(angle) sin(angle) 0];
            obj.direction = Vector3D;
            obj.direction.setValue(dir);
        end
        
        function computeWeakStiffTensor(obj)
            E1  = 1;
            nu1 = 0;
            E0  = 1e-1*E1;
            nu0 = 0;
            obj.stiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.weakTensor  = IsotropicConstitutiveTensor(E0,nu0);
        end
        
        function computeHomogVoigtPlaneStressTensor(obj)
            c0     = obj.weakTensor;
            c1     = obj.stiffTensor;
            dir{1} = obj.direction;
            m1     = 1;
            seqHomog = HomogVoigtPlaneStressHomogenizer(c0,c1,dir,m1,obj.theta);
            obj.hvpTensor  = seqHomog.getPlaneStressHomogenizedTensor();
        end
        
        function computeVoigtHomogPlaneStressTensor(obj)
            c0     = obj.weakTensor;
            c1     = obj.stiffTensor;
            dir{1} = obj.direction;
            m1     = 1;
            seqHomog = VoigtHomogPlaneStressHomogenizer(c0,c1,dir,m1,obj.theta);
            obj.vhpTensor  = seqHomog.getPlaneStressHomogenizedTensor();
        end
        
    end

end

