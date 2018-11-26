classdef testAmplificatorTensorNumericVsExplicitForSeqLam < testShowingError
    
    
    properties (Access = private)
        ampTensorNum
        ampTensorExp
        fiberDirection
    end
    
    properties (Access = protected)
       tol = 1e-12;
    end
       
    
    methods (Access = public)
        
        function obj = testAmplificatorTensorNumericVsExplicitForSeqLam() 
            obj.createNumericalAmplificationTensor()
            obj.createExplicitAmplificationTensor()
        end
      
        
    end
    
    methods (Access = private)
        
        function createFiberDirection(obj)
           dir  = [1 0 0];
           obj.fiberDirection = Vector3D;
           obj.fiberDirection.setValue(dir);
           obj.fiberDirection.normalize();
        end
                                
        function createNumericalAmplificationTensor(obj)
            obj.createFiberDirection()
            dir            = obj.fiberDirection;
            LevelOfFibers  = 3;
            FamilyName = 'HorizontalLaminate';
            LevelStr   = num2str(LevelOfFibers);
            name = strcat(FamilyName,LevelStr);
            printTopology  = false;
            iter           = 0;
            homogenizer    = NumericalFiberHomogenizer(dir,...
                             LevelOfFibers,name,...
                             printTopology,iter);
            Ch             = homogenizer.getCh();
            Volume         = homogenizer.getVolume();
            P              = homogenizer.getAmplificatorTensor();
            obj.ampTensorNum = SymmetricFourthOrderPlaneStressVoigtTensor();
            obj.ampTensorNum.createRandomTensor();
        end
        
        function createExplicitAmplificationTensor(obj)
            obj.ampTensorExp = SymmetricFourthOrderPlaneStressVoigtTensor();
            obj.ampTensorExp.createRandomTensor();
        end

    end    
    
    methods (Access = protected)
        
        function computeError(obj)
            va = obj.ampTensorNum.getValue();
            ta = obj.ampTensorExp.getValue();
            obj.error = norm(ta(:) - va(:))/norm(ta(:));
        end
        
    end
    

end

