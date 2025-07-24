classdef Homogenizer < handle
    
    properties (Access = private)
        anisotropicTensors
        laminateParams
        theta
        rank
        
        matTensor
        fibTensor
        homogTensor
        aniTensor
        incremTensor
        
        invMatTensor
        invFibTensor
        invThetaTensor
    end
    
    methods (Access = public)
        
        function obj = Homogenizer(C0,C1,Ani,mi,theta)
            obj.init(C0,C1,Ani,mi,theta)
            obj.computeIncrementalTensor();
            obj.computeAnisotropicContributionTensor()
            obj.computeThetaTensor();
            obj.computeHomogenizedTensor();
        end
        
        function t = getHomogenizedTensor(obj)
            t = obj.homogTensor;
        end
         
    end
    
    methods (Access = private)
        
        function init(obj,C0,C1,Ani,mi,theta)
            obj.matTensor = C0;
            obj.fibTensor = C1;
            obj.anisotropicTensors = Ani;
            obj.laminateParams     = mi;
            obj.theta              = theta;
            obj.rank               = length(obj.laminateParams);
            obj.computeMatrixAndFiberComplianceTensors();
        end
        
        function computeMatrixAndFiberComplianceTensors(obj)
            obj.invMatTensor = Inverter.invert(obj.matTensor);
            obj.invFibTensor = Inverter.invert(obj.fibTensor);
        end
        
        function computeIncrementalTensor(obj)
            S0  = obj.invMatTensor;
            S1  = obj.invFibTensor;
            S01 = S0.clone();
            s0 = S0.getValue();
            s1 = S1.getValue();
            s01 = s0 - s1;
            S01.setValue(s01);
            C01 = Inverter.invert(S01);
            obj.incremTensor = C01;
        end
        
        function computeAnisotropicContributionTensor(obj)
            tensorSize = obj.anisotropicTensors{1}.getTensorSize();
            Ca = zeros(tensorSize);
            for irank = 1:obj.rank
                Cai = obj.anisotropicTensors{irank}.getValue();
                mi  = obj.laminateParams(irank);
                Ca  = Ca + mi*Cai;
            end
            obj.aniTensor = obj.anisotropicTensors{1}.clone();
            obj.aniTensor.setValue(Ca);
        end
        
        function computeThetaTensor(obj)
            C01 = obj.incremTensor;
            Ca  = obj.aniTensor;
            Ctheta = C01.clone();
            c01 = C01.getValue();
            ca  = Ca.getValue();
            ctheta = c01 + obj.theta*ca;
            Ctheta.setValue(ctheta);
            Stheta = Inverter.invert(Ctheta);
            obj.invThetaTensor = Stheta;
        end
        
        function computeHomogenizedTensor(obj)
            S1     = obj.invFibTensor;
            Stheta = obj.invThetaTensor;
            Sh = S1.clone();
            s1 = S1.getValue();
            stheta = Stheta.getValue();
            sh = s1 + (1-obj.theta)*stheta;
            Sh.setValue(sh);
            Ch = Inverter.invert(Sh);
            obj.homogTensor = Ch;
        end
        
    end
    
end

