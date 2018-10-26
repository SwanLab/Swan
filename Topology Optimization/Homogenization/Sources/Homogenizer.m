classdef Homogenizer < handle
    
    properties (Access = private)
        anisotropicTensors
        laminateParams
        Theta
        Rank
        
        MatTensor
        FibTensor
        HomogTensor
        AniTensor
        IncremTensor
        
        InvMatTensor
        InvFibTensor
        InvThetaTensor
    end
    
    methods (Access = public)
        
        function obj = Homogenizer(C0,C1,Ani,mi,theta)
            obj.init(C0,C1,Ani,mi,theta)
            obj.computeIncrementalTensor();
            obj.computeAnisotropicContributionTensor()
            obj.computeThetaTensor();
            obj.computeHomogenizedTensor();
        end
        
        function T = getHomogenizedTensor(obj)
            T = obj.HomogTensor;
        end
         
    end
    
    methods (Access = private)
        
        function init(obj,C0,C1,Ani,mi,theta)
            obj.MatTensor = C0;
            obj.FibTensor = C1;
            obj.anisotropicTensors = Ani;
            obj.laminateParams     = mi;
            obj.Theta              = theta;
            obj.Rank               = length(obj.laminateParams);
            obj.computeMatrixAndFiberComplianceTensors();
        end
        
        function computeIncrementalTensor(obj)
            S0  = obj.InvMatTensor;
            S1  = obj.InvFibTensor;
            S01 = S0 - S1;
            C01 = obj.invertTensor(S01);
            obj.IncremTensor = C01;
        end
            
        function computeMatrixAndFiberComplianceTensors(obj)
            obj.InvMatTensor = obj.invertTensor(obj.MatTensor.getValue());
            obj.InvFibTensor = obj.invertTensor(obj.FibTensor.getValue());
        end
        
        function computeAnisotropicContributionTensor(obj)
            TensorSize = size(obj.anisotropicTensors{1}.getValue());
            Ca = zeros(TensorSize);
            for irank = 1:obj.Rank
                Cai = obj.anisotropicTensors{irank}.getValue();
                mi  = obj.laminateParams(irank);
                Ca  = Ca + mi*Cai;
            end
            obj.AniTensor = Ca;
        end
        
        function computeThetaTensor(obj)
            C01 = obj.IncremTensor;
            Ca  = obj.AniTensor;
            Ctheta = C01 + obj.Theta*Ca;
            Stheta = obj.invertTensor(Ctheta);
            obj.InvThetaTensor = Stheta;
        end
        
        function computeHomogenizedTensor(obj)
            S1     = obj.InvFibTensor;
            Stheta = obj.InvThetaTensor;
            Sh = S1 + (1-obj.Theta)*Stheta;
            Ch = obj.invertTensor(Sh);
            obj.HomogTensor = Ch;
        end
    end
    
    methods (Access = private, Static)
        
        function invT = invertTensor(TensorValue)
            T = FourthOrderTensor();
            T.setValue(TensorValue);
            invT = Inverter.invert(T);
        end
        
    end
end

