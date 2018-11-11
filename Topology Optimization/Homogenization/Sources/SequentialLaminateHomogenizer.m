classdef SequentialLaminateHomogenizer < handle
    
    
    properties (Access = public)
       HomogenizedTensor
    end
    
    properties (Access = private)
        Directions
        LaminateParams
        Rank
        Theta
        FiberTensor
        MatrixTensor
        AnisotropicTensors
        IncrementalTensor
        AniTensor
        ThetaTensor   
    end
    
    methods (Access = public)
        
        function obj = SequentialLaminateHomogenizer(FiberTensor,MatrixTensor,directions,LaminateParams,Theta)
            obj.init(FiberTensor,MatrixTensor,directions,LaminateParams,Theta);
            obj.computeAnisotropicTensors();
            obj.computeIncrementalTensor();
            obj.computeAnisotropicContributionTensor()
            obj.computeThetaTensor();
            obj.computeHomogenizedTensor();
            obj.computePlaneStressHomogenizedTensor();
        end
    end
    
    methods (Access = private)
        function init(obj,FiberTensor,MatrixTensor,directions,LaminateParams,Theta)
            obj.FiberTensor    = FiberTensor;
            obj.MatrixTensor   = MatrixTensor;
            obj.LaminateParams = LaminateParams;
            obj.Rank           = length(obj.LaminateParams);
            obj.Theta          = Theta;
            obj.normalizeDirections(directions)
            obj.IncrementalTensor = fourthOrderTensor();
            obj.AniTensor         = fourthOrderTensor();
            obj.ThetaTensor       = fourthOrderTensor();
            obj.HomogenizedTensor = fourthOrderTensor();
            obj.computeMatrixAndFiberComplianceTensors();
            obj.checkConsistencyLaminateParamsAndDirections()
        end
        
        function normalizeDirections(obj,directions)
            for irank = 1:obj.Rank
               DirectionNorm =  norm(directions(irank,:));
               obj.Directions(irank,:) = directions(irank,:)/DirectionNorm;
            end
        end
        
        function computeAnisotropicTensors(obj)
            for irank = 1:obj.Rank
            obj.AnisotropicTensors{irank} = AnisotropicContributionTensor(obj.FiberTensor,obj.Directions(irank,:));
            obj.AnisotropicTensors{irank}.computeTensorVoigtInPlaneStress()     
            end
        end
        
        function computeMatrixAndFiberComplianceTensors(obj)
            obj.MatrixTensor.InverseTensorVoigt = inv(obj.MatrixTensor.tensorVoigt);
            obj.FiberTensor.InverseTensorVoigt  = inv(obj.FiberTensor.tensorVoigt);         
        end
        
        function checkConsistencyLaminateParamsAndDirections(obj)
            assert(size(obj.Directions,1) == length(obj.LaminateParams))
        end
        
        function computeIncrementalTensor(obj)
            S0  = obj.MatrixTensor.InverseTensorVoigt;
            S1  = obj.FiberTensor.InverseTensorVoigt;
            S01 = S0 - S1;
            C01 = inv(S01);
            obj.IncrementalTensor.InverseTensorVoigt =  S01;
            obj.IncrementalTensor.tensorVoigt = C01;           
        end
        
        function computeAnisotropicContributionTensor(obj)
            VoigtSize = size(obj.AnisotropicTensors{1}.tensorVoigt); 
            Ca = zeros(VoigtSize);
            for irank = 1:obj.Rank
                 CaThisRank       = obj.AnisotropicTensors{irank}.tensorVoigt;
                 LamParamThisRank = obj.LaminateParams(irank);
                 Ca = Ca + LamParamThisRank*CaThisRank;
            end                        
            obj.AniTensor.tensorVoigt = Ca;
        end
        
        function computeThetaTensor(obj)
            C01 = obj.IncrementalTensor.tensorVoigt;
            Ca  = obj.AniTensor.tensorVoigt;
            Ctheta = C01 + obj.Theta*Ca;
            Stheta = inv(Ctheta);
            obj.ThetaTensor.tensorVoigt = Ctheta ;
            obj.ThetaTensor.InverseTensorVoigt = Stheta;
        end
        
        function computeHomogenizedTensor(obj)
            S1 = obj.FiberTensor.InverseTensorVoigt;
            Stheta = obj.ThetaTensor.InverseTensorVoigt;
            Sh = S1 + (1-obj.Theta)*Stheta;
            Ch = inv(Sh);
            obj.HomogenizedTensor.InverseTensorVoigt = Sh ;
            obj.HomogenizedTensor.tensorVoigt = Ch;
        end
        
        function computePlaneStressHomogenizedTensor(obj)
            obj.HomogenizedTensor.computeTensorVoigtInPlaneStress();
        end

    end

    


      
end

