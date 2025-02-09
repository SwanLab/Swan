classdef testCommutingHomogPlaneStress < handle
    
    properties (Access = private)
        theta
        direction
        stiffTensor
        weakTensor
    end
    
    properties (Access = protected)
        vhpTensor
        vphTensor
    end
    
    methods (Access = public)
        
        function obj = testCommutingHomogPlaneStress
            obj.init()
            obj.computeWeakStiffTensor()
            obj.computeVoigtPlaneStressHomog()
            obj.computeVoigtHomogPlaneStressTensor()
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
            nu1 = obj.createPoissonValue;
            E0  = 1e-3*E1;
            nu0 = nu1;
            obj.stiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.weakTensor  = IsotropicConstitutiveTensor(E0,nu0);
        end
        
        function computeVoigtPlaneStressHomog(obj)
            c0     = obj.weakTensor;
            c1     = obj.stiffTensor;
            dir{1} = obj.direction;
            m1     = 1;
            seqHomog = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,m1,obj.theta);
            obj.vphTensor  = seqHomog.getPlaneStressHomogenizedTensor();
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

    methods (Access = protected, Abstract, Static)
        createPoissonValue(obj)
    end
    
end