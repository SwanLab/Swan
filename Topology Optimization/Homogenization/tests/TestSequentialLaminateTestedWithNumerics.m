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
            C0 = obj.WeakTensor;
            C1 = obj.StiffTensor;
            dir = obj.LaminateDirection;
            m1 = 1;
            SeqHomog = VoigtHomogPlaneStressHomogenizer(C0,C1,dir,m1,obj.Theta);  
            obj.SeqLamCh  = SeqHomog.getPlaneStressHomogenizedTensor();
      end
        
        function computeRank2HomogenizerTensor(obj)
            
            C0 = obj.WeakTensor;
            C1 = obj.StiffTensor;
            dir = obj.LaminateDirection;
            m1 = 1;
            SeqHomog = VoigtPlaneStressHomogHomogenizer(C0,C1,dir,m1,obj.Theta);  
            obj.Rank2Ch  = SeqHomog.getPlaneStressHomogenizedTensor();                        
        end
                 
        function computeMixtureTheoryTensor(obj)
           C1 = obj.StiffTensor;
           C0 = obj.WeakTensor; 
           Dir = [0 0 1];
           Angle = -acos(dot(obj.FiberDirection,[1 0 0]));
           Vfrac = obj.Theta;
           Homogenizer = MixtureTheoryHomogenizer(C1,C0,Dir,Angle,Vfrac);
           obj.MixtureCh = Homogenizer.Ch;            
        end
  
 
        function hasPassed = hasPassed(obj)
            ChSL    = double(obj.SeqLamCh);
            ChNum   = double(obj.NumericalCh); 
            ChMix   = double(obj.MixtureCh); 
            ChRank  = double(obj.Rank2Ch);  
            
            firstCondition  = obj.relativeNorm(ChSL,ChNum)   < 1e-2;
            secondCondition = obj.relativeNorm(ChMix,ChNum)  < 1e-3;
            thirdCondition  = obj.relativeNorm(ChRank,ChNum) < 1e-10;
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

