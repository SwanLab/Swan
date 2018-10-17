classdef TestSequentialLaminateTestedWithNumerics < test

  
    properties (Access = protected)
        LaminateDirection
        Theta
        FractionVolume
       
        StiffTensor 
        WeakTensor
        MaterialValues
       
        MixtureCh
        Rank2Ch
        NumericalCh
        SeqLamCh          
        
        FiberDirection
    end
    
  
    methods (Access = protected)

        function compute(obj)
            obj.init()
            obj.computeNumericallyChForLaminate();
            obj.loadFractionVolume()
            obj.computeWeakAndStiffTensorsFromNumericalHomogenizerData();
            obj.computeSequentialLaminateTensor();
            obj.computeRank2HomogenizerTensor();
            obj.computeMixtureTheoryTensor();
        end
        
        function init(obj)
           obj.loadLaminateDirection()
           obj.loadFiberDirection()
        end
        
        function computeNumericallyChForLaminate(obj)   
            OutPutNameFile     = 'SeqLaminate';
            LevelOfFibers      = 3;
            PrintTopology      = true;
            NumHomogenizer     = NumericalFiberHomogenizer(...
                                 obj.FiberDirection,LevelOfFibers,...
                                 OutPutNameFile,PrintTopology);
            obj.NumericalCh    = NumHomogenizer.Ch;
            obj.MaterialValues = NumHomogenizer.MaterialValues;
            obj.FractionVolume = NumHomogenizer.Volume;
        end

        function loadFractionVolume(obj)
           obj.Theta = obj.FractionVolume;
        end
        
        function computeWeakAndStiffTensorsFromNumericalHomogenizerData(obj)
            E1  = obj.MaterialValues.E_plus;
            nu1 = obj.MaterialValues.nu_plus;
            E0  = obj.MaterialValues.E_minus;
            nu0 = obj.MaterialValues.nu_minus;            
            obj.StiffTensor = obj.createIsotropicTensor(E1,nu1);
            obj.WeakTensor  = obj.createIsotropicTensor(E0,nu0);
        end
        
        function Tensor = createIsotropicTensor(obj,E,nu)
            Tensor = IsotropicConstitutiveTensor3D(E,nu);
            Tensor.computeTensorVoigtInPlaneStress()
        end
        
        function computeSequentialLaminateTensor(obj)
            C1 = obj.StiffTensor;
            C0 = obj.WeakTensor;
            dir = obj.LaminateDirection;
            m1 = 1;
            SeqHomogenizer = SequentialLaminateHomogenizer(C1,C0,dir,m1,obj.Theta);                                              
            Ch  = SeqHomogenizer.HomogenizedTensor;
            obj.SeqLamCh = double(Ch.tensorVoigtInPlaneStress);          
       end
        
        function computeRank2HomogenizerTensor(obj)
            Params = [1 0];
            Dir = [obj.LaminateDirection;obj.LaminateDirection];
            C1 = obj.StiffTensor;
            C0 = obj.WeakTensor;
            Rank2 = RankTwoLaminateHomogenizer(C1,C0,Dir,Params,obj.Theta);
            obj.Rank2Ch = double(Rank2.Ch);
        end
                 
        function computeMixtureTheoryTensor(obj)
           C1 = obj.StiffTensor;
           C0 = obj.WeakTensor; 
           Homogenizer = MixtureTheoryHomogenizer(C1,C0,obj.Theta);
           obj.MixtureCh = Homogenizer.Ch;            
        end
  
 
        function hasPassed = hasPassed(obj)
            ChSL    = double(obj.SeqLamCh);
            ChNum   = double(obj.NumericalCh); 
            ChMix   = double(obj.MixtureCh); 
            ChRank2 = double(obj.Rank2Ch);            
            firstCondition  = obj.relativeNorm(ChNum,ChSL)   < 1e-2;
            secondCondition = obj.relativeNorm(ChMix,ChSL)   < 1e-2;
            thirdCondition  = obj.relativeNorm(ChRank2,ChSL) < 3*1e-2;
            hasPassed = firstCondition & secondCondition & thirdCondition;
        end
        
        function RelNorm = relativeNorm(obj,A,B)
            RelNorm = norm(A - B)/norm(B);
        end
        
    end
    
    methods (Abstract,Access = protected)
        loadLaminateDirection(obj)
        loadFiberDirection(obj)
    end
end

