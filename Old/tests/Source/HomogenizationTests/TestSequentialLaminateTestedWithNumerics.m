classdef TestSequentialLaminateTestedWithNumerics < handle

    properties (Access = private)
        microFile = 'RVE_Square_Triangle_FineFine';
        fileOutputName = 'SeqLaminate';
    end

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

    methods (Access = public)

        function hasPassed = hasPassed(obj) 
            ChNum   = obj.NumericalCh.getValue();
            ChSL    = obj.SeqLamCh.getValue();
            ChMix   = obj.MixtureCh.getValue();
            ChRank  = obj.Rank2Ch.getValue();
            firstCondition  = obj.relativeNorm(ChSL,ChNum)   < 1e-2;
            secondCondition = obj.relativeNorm(ChMix,ChNum)  < 1e-3;
            thirdCondition  = obj.relativeNorm(ChRank,ChNum) < 1e-10;
            hasPassed = firstCondition & secondCondition & thirdCondition;
        end

    end

    methods (Access = protected)

        function compute(obj)
            obj.init();
            obj.computeNumericallyChForLaminate();
            obj.loadFractionVolume()
            obj.computeWeakAndStiffTensorsFromNumericalHomogenizerData();
            obj.computeSequentialLaminateTensor();
            obj.computeRank2HomogenizerTensor();
            obj.computeMixtureTheoryTensor();
        end

    end
    
    methods (Access = private)
                
        function init(obj)
           obj.loadLaminateDirection()
           obj.loadFiberDirection()
        end
        
        function computeNumericallyChForLaminate(obj)
           d = obj.createNumericalHomogenizerDataBase();
           homog = NumericalHomogenizer(d);
           homog.compute();
           obj.NumericalCh    = obj.rotateCh(homog);
           obj.MaterialValues = homog.matValues;
           obj.FractionVolume = homog.cellVariables.volume;
        end
        
        function d = createNumericalHomogenizerDataBase(obj)
            nDB = NumericalHomogenizerDataBase(obj.microFile);
            d = nDB.dataBase;
            d.outFileName = obj.fileOutputName;
            d.hasToCaptureImage = false;
        end
        
        function loadFractionVolume(obj)
           obj.Theta = obj.FractionVolume;
        end
        
        function computeWeakAndStiffTensorsFromNumericalHomogenizerData(obj)
            E1  = obj.MaterialValues.E_plus;
            nu1 = obj.MaterialValues.nu_plus;
            E0  = obj.MaterialValues.E_minus;
            nu0 = obj.MaterialValues.nu_minus;
            obj.StiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.WeakTensor  = IsotropicConstitutiveTensor(E0,nu0);
        end
        
        function computeSequentialLaminateTensor(obj)
            C0 = obj.WeakTensor;
            C1 = obj.StiffTensor;
            dir{1} = obj.LaminateDirection;
            m1 = 1;
            SeqHomog = VoigtHomogPlaneStressHomogenizer(C0,C1,dir,m1,obj.Theta);
            obj.SeqLamCh  = SeqHomog.getPlaneStressHomogenizedTensor();
        end
        
        function computeRank2HomogenizerTensor(obj)
            C0 = obj.WeakTensor;
            C1 = obj.StiffTensor;
            dir{1} = obj.LaminateDirection;
            m1 = 1;
            SeqHomog = VoigtPlaneStressHomogHomogenizer(C0,C1,dir,m1,obj.Theta);
            obj.Rank2Ch  = SeqHomog.getPlaneStressHomogenizedTensor();
        end
        
        function computeMixtureTheoryTensor(obj)
            C1 = obj.StiffTensor;
            C0 = obj.WeakTensor;
            d = [0 0 1];
            dir = Vector3D;
            dir.setValue(d);
            dir.normalize();
            angle = -acos(dot(obj.FiberDirection.getValue(),[1 0 0]));
            vFrac = obj.Theta;
            homogenizer = MixtureTheoryHomogenizer(C1,C0,dir,angle,vFrac);
            obj.MixtureCh = homogenizer.Ch;
        end
        
        function Ch = rotateCh(obj,homog)
            dir = obj.FiberDirection;
            r = ChRotatorForFiberHomogenizer();
            Ch = r.rotate(dir,homog.cellVariables.Ch());
        end
       
    end
    
    methods (Access = private, Static)
        
        function relNorm = relativeNorm(A,B)
            relNorm = norm(A - B)/norm(B);
        end

    end
   
    methods (Abstract,Access = protected)
        loadLaminateDirection(obj)
        loadFiberDirection(obj)
    end
end