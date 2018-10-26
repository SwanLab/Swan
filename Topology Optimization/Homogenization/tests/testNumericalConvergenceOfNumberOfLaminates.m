classdef testNumericalConvergenceOfNumberOfLaminates < test
    
    properties (Access = private)
        
        FiberDirection
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
           obj.FiberDirection        = [1 0 0];
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
            Dir            = obj.FiberDirection;
            PrintTopology  = true;
            Homogenizer    = NumericalFiberHomogenizer(Dir,LoF,obj.Name,...
                             PrintTopology);
            obj.Ch         = Homogenizer.Ch;  
            obj.Volume     = Homogenizer.Volume;
        end

        function computeChNorm(obj)
            for iCh = 1:size(obj.Ch,2)
                obj.ChNorm(iCh) = norm(obj.AllCh{iCh});
            end
        end
        
        function AreChSimilar = ComputeChSimilarity(obj)
            NormCh = obj.ChNorm;
            meanChNorm = mean(NormCh); 
            AreChSimilar = norm(NormCh - meanChNorm) < 1e-12;
        end
        
        function AreVolumesSimilar = ComputeVolumeSimilarity(obj)
            Volumes = obj.AllVolume;
            meanVolumes = mean(Volumes);
            AreVolumesSimilar = norm(Volumes - meanVolumes) < 1e-12;
        end
        
    end
    
    methods (Access = protected)
        
        function TestHasPassed = hasPassed(obj)
            AreChSimiliar    = obj.ComputeChSimilarity();
            AreVolumeSimilar = obj.ComputeVolumeSimilarity();
            TestHasPassed = AreChSimiliar && AreVolumeSimilar;
        end
        
    end
end

