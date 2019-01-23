classdef testHomogenizationLaminateForAnisotropicTensors < ...
         testShowingError
    
     
     properties (Access = protected)
         testName = 'HomogenizationLaminateForAnisotropicTensors';
         tol = 1e-6;
     end
     
     properties (Access = private)
         stiffTensor
         weakTensor
         ChForIso
         ChForAni
         laminateDirection
         theta
     end
     
     methods (Access = public)
         
         function obj = testHomogenizationLaminateForAnisotropicTensors()
             obj.init()
             obj.computeHomogenizerForIsotropicMaterials()
             obj.computeHomogenizerForAnisotropicMaterials()
         end
         
     end
     
     
     methods (Access = protected)
         
         function init(obj)
             obj.createTensors()
             obj.createLaminateDirection()
             obj.theta = 0.8;
         end
         
        function createLaminateDirection(obj)
            d = [0 1 0];            
            dir = Vector3D;
            dir.setValue(d);
            dir.normalize()
            obj.laminateDirection = dir;
        end
         
         function createTensors(obj)
            E1 = 1;
            E0 = 1e-3;
            nu1 = 1/3;
            nu0 = 1/3;
            obj.stiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.weakTensor  = IsotropicConstitutiveTensor(E0,nu0);
         end
         
         function computeHomogenizerForIsotropicMaterials(obj)
            C0 = obj.weakTensor;
            C1 = obj.stiffTensor;
            dir{1} = obj.laminateDirection;
            m1 = 1;
            SeqHomog = VoigtPlaneStressHomogHomogenizer(C0,C1,dir,m1,obj.theta);
            obj.ChForIso  = SeqHomog.getPlaneStressHomogenizedTensor();
         end
         
         function computeHomogenizerForAnisotropicMaterials(obj)
            C0 = obj.weakTensor;
            C1 = obj.stiffTensor;
            dir = obj.laminateDirection;
            Lam = AnisotropicLaminateHomogenizer(C0,C1,dir,obj.theta); 
            obj.ChForAni = Lam.getHomogenizedTensor();
         end
         
         
         function computeError(obj)
             Ch4Iso = obj.ChForIso.getValue();
             Ch4Ani = obj.ChForAni.getValue();
             obj.error = norm(Ch4Iso - Ch4Ani);
         end
         
     end
     
end

