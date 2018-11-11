classdef TestGeneralTwoRankSequentialLaminate < test
    
    
    
    properties (Access = private)

      FractionVolume      
      
      Directions      
      LamParams
      
      StiffTensor
      WeakTensor
      
      Rank2Ch
      SeqLamCh
      
    end
    
    methods (Access = public) 
        
        function obj = TestGeneralTwoRankSequentialLaminate()
            obj.init()
            obj.computeTwoRankSequentialLaminate()
            obj.computeGeneralTwoRankSequentialLaminate()
        end

    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            RankTwoCh  = double(obj.Rank2Ch);
            SqCh = double(obj.SeqLamCh);
            hasPassed = norm(RankTwoCh - SqCh)/norm(SqCh) < 1e-2;
        end
    end
    
    methods (Access = private)
        
        function init(obj)            
            obj.FractionVolume = 0.8000;
            obj.LamParams = [1 0];
            obj.loadLaminateDirections()
            obj.createStiffAndWeakTensors()
        end
              
        function loadLaminateDirections(obj)
            obj.Directions(1,:) = [1     0     0];
            obj.Directions(2,:) = [1     3     0];
            obj.Directions = obj.normalizeDirections(obj.Directions);
        end       
        
        function createStiffAndWeakTensors(obj)
            epsilon = 1.0000e-03;
            E1  = 1;
            nu1 = 1/3;
            E0  = epsilon*E1;
            nu0 = 1/3;           
            obj.StiffTensor = obj.createIsotropicTensor(E1,nu1);
            obj.WeakTensor  = obj.createIsotropicTensor(E0,nu0);
        end
        
        function Tensor = createIsotropicTensor(obj,E,nu)
            Tensor = IsotropicConstitutiveTensor3D(E,nu);
            Tensor.computeTensorVoigtInPlaneStress()
        end
            
        function computeTwoRankSequentialLaminate(obj)           
            Param    = obj.LamParams;
            Dir      = obj.Directions;
            Theta = obj.FractionVolume;
            C1    = obj.StiffTensor;
            C0    = obj.WeakTensor;
            
            Homogenizer = RankTwoLaminateHomogenizer(C1,C0,Dir,Param,Theta);
            obj.Rank2Ch = Homogenizer.Ch;
        end
        
        function computeGeneralTwoRankSequentialLaminate(obj)            
            C1       = obj.StiffTensor;
            C0       = obj.WeakTensor;
            Param    = obj.LamParams;
            Dir      = obj.Directions;
            Theta    = obj.FractionVolume;            
                
            SeqHomog     = SequentialLaminateHomogenizer(C1,C0,Dir,Param,Theta);
            HomogTensor  = SeqHomog.HomogenizedTensor;
            obj.SeqLamCh = HomogTensor.tensorVoigtInPlaneStress;
  
        end

    end
    
    methods (Access = private, Static)
       
        function NormalizedDir = normalizeDirections(directions)
            NormalizedDir = zeros(size(directions));
            rank = size(directions,1);
            for irank = 1:rank
                DirectionNorm =  norm(directions(irank,:));
                NormalizedDir(irank,:) = directions(irank,:)/DirectionNorm;
            end
        end
        
    end
    
    
    
end