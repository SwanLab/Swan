classdef testNumericalConvergenceOfNumberOfLaminates < handle
    
    properties (Access = public)
       tol = 1e-12 
    end
    
    properties (Access = private)
        fiberDirection
        nLevelsOfFibers
        Ch
        outFileName
        
        AllCh
        AllVolume
        
        Volume
        ChNorm
    end
    
    methods (Access = public)
        
        function obj = testNumericalConvergenceOfNumberOfLaminates()
            obj.createFiberDirection()
            obj.createNumberOfLevelsOfFibers()
            obj.computeAllChAndVolumes()
            obj.computeChNorm()
        end
        
        function error = computeError(obj)
            errorCh = obj.ComputeChSimilarity();
            errorVol = obj.ComputeVolumeSimilarity();
            error = max(errorCh,errorVol);
        end
        
    end
    
    methods (Access = private)
        function createFiberDirection(obj)
           dir  = [1 0 0];
           obj.fiberDirection = Vector3D;
           obj.fiberDirection.setValue(dir);
           obj.fiberDirection.normalize();
        end
        
        function createNumberOfLevelsOfFibers(obj)
            obj.nLevelsOfFibers = 3;
        end
        
        function computeAllChAndVolumes(obj)
            nlevel = obj.nLevelsOfFibers;
            for iLevel = 1:nlevel
                obj.computeName(iLevel);
                obj.computeHomogenization(iLevel);
                obj.AllCh{iLevel}      = obj.Ch;
                obj.AllVolume(iLevel,:)  = obj.Volume;
            end
        end
        
        function computeName(obj,iLevel)
            familyName = 'HorizontalLaminate';
            LevelStr   = num2str(iLevel);
            obj.outFileName = strcat(familyName,LevelStr);
        end
        
        function computeHomogenization(obj,LoF)
           d = obj.computeNumericalHomogenizerDataBase(LoF);
           homog = NumericalHomogenizer(d);
           homog.compute();
           obj.Ch = obj.rotateCh(homog);
           obj.Volume = homog.cellVariables.volume;
        end
        
        function d = computeNumericalHomogenizerDataBase(obj,LoF)
            microFile = 'RVE_Square_Triangle_FineFine';
            nDB = NumericalHomogenizerDataBase(microFile);
            d = nDB.dataBase;
            d.outFileName = obj.outFileName;
            d.levelSetDataBase.levFib = LoF;
            d.hasToCaptureImage = false;
        end
        
        function C = rotateCh(obj,homog)
           dir   = obj.fiberDirection;
           r = ChRotatorForFiberHomogenizer();
           C = r.rotate(dir,homog.cellVariables.Ch());
        end

        function computeChNorm(obj)
            for iCh = 1:numel(obj.AllCh)
                obj.ChNorm(iCh) = norm(obj.AllCh{iCh}.getValue());
            end
        end
        
        function ChError = ComputeChSimilarity(obj)
            NormCh = obj.ChNorm;
            meanChNorm = mean(NormCh);
            ChError = norm(NormCh - meanChNorm);
        end
        
        function volError = ComputeVolumeSimilarity(obj)
            Volumes = obj.AllVolume;
            meanVolumes = mean(Volumes);
            volError = norm(Volumes - meanVolumes);
        end
        
    end

end

