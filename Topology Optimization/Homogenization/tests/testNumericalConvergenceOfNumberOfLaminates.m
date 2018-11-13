classdef testNumericalConvergenceOfNumberOfLaminates < testShowingError
    
    properties (Access = protected)
       tol = 1e-12 
    end
    
    properties (Access = private)
        
        fiberDirection
        NumberOfLevelsOfFibers
        Ch
        Name
        
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
        
    end
    
    methods (Access = private)
        function createFiberDirection(obj)
           dir  = [1 0 0];
           obj.fiberDirection = Vector3D;
           obj.fiberDirection.setValue(dir);
           obj.fiberDirection.normalize();
        end
        
        function createNumberOfLevelsOfFibers(obj)
            obj.NumberOfLevelsOfFibers = 4;
        end
        
        function computeAllChAndVolumes(obj)
            nlevel = obj.NumberOfLevelsOfFibers;
            for iLevel = 1:nlevel
                obj.computeName(iLevel);
                obj.computeHomogenization(iLevel);
                obj.AllCh{iLevel}      = obj.Ch;
                obj.AllVolume(iLevel)  = obj.Volume;
            end
        end
        
        function computeName(obj,iLevel)
            FamilyName = 'HorizontalLaminate';
            LevelStr   = num2str(iLevel);
            obj.Name = strcat(FamilyName,LevelStr);
        end
        
        function computeHomogenization(obj,LoF)
            Dir            = obj.fiberDirection;
            PrintTopology  = true;
            Homogenizer    = NumericalFiberHomogenizer(Dir,LoF,obj.Name,...
                             PrintTopology);
            obj.Ch         = Homogenizer.Ch;  
            obj.Volume     = Homogenizer.Volume;
        end

        function computeChNorm(obj)
            for iCh = 1:size(obj.Ch,2)
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
    
    methods (Access = protected)
        
        function computeError(obj)
            errorCh = obj.ComputeChSimilarity();
            errorVol = obj.ComputeVolumeSimilarity();
            obj.error = max(errorCh,errorVol);
        end
        
    end
end

