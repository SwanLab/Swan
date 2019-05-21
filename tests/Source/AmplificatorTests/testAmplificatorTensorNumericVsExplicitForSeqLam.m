classdef testAmplificatorTensorNumericVsExplicitForSeqLam < testShowingError
    
    
    properties (Access = private)
        ampTensorNum
        ampTensorExp
        fiberDirection
        matValues
        C1
        C0
        Ch
        theta
    end
    
    properties (Access = protected)
       tol = 1e-10;
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
            LevelOfFibers  = 4;
            FamilyName = 'HorizontalLaminate';
            LevelStr   = num2str(LevelOfFibers);
            name = strcat(FamilyName,LevelStr);
            printTopology  = true;
            iter           = 0;
            homogenizer    = NumericalFiberHomogenizer(dir,...
                             LevelOfFibers,name,...
                             printTopology,iter);
            obj.Ch         = homogenizer.getCh();
            obj.theta      = homogenizer.getVolume();
            P              = homogenizer.getAmplificatorTensor();
            obj.matValues  = homogenizer.getMaterialValues();
            obj.ampTensorNum = P;
        end
        
        function createExplicitAmplificationTensor(obj)
            obj.createMaterialTensors()
            c0 = obj.C0;
            c1 = obj.C1;
            dir = obj.fiberDirection;
            dirVal = dir.getValue;
            dirVal = [-dirVal(2),dirVal(1),0];
            dir.setValue(dirVal);
            t = obj.theta;
            Lam = AnisotropicLaminateHomogenizer(c0,c1,dir,t);
            P = Lam.getAmplificatorTensor();
            Ptens = CompliancePlaneStressTensor();
            Ptens.setValue(P);
            PtensV = Tensor2VoigtConverter.convert(Ptens);
            obj.ampTensorExp = PtensV;
        end
        
        function createMaterialTensors(obj)
            E1  = obj.matValues.E_plus;
            nu1 = obj.matValues.nu_plus;
            E0  = obj.matValues.E_minus;
            nu0 = obj.matValues.nu_minus;
            obj.C1  = IsotropicConstitutiveTensor(E1,nu1);
            obj.C0  = IsotropicConstitutiveTensor(E0,nu0);
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

