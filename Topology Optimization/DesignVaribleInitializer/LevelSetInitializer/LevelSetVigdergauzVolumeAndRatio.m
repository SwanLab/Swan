classdef LevelSetVigdergauzVolumeAndRatio < LevelSetCreator
    
    properties (Access = private)
        ratio
        volume
    end
    
    methods (Access = public)
        
        function obj = LevelSetVigdergauzVolumeAndRatio(cParams)
            obj.init(cParams);
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            volum = obj.volume;
            phi = atan(obj.ratio);
            if obj.isMicroStructureValid(phi,volum)
                cParams.x = obj.nodeCoord(:,1);                
                cParams.y = obj.nodeCoord(:,2);  
                cParams.volum = volum;
                cParams.phi = phi;
                lsVig = LevelSetLipung(cParams);
                obj.levelSet = lsVig.value;
            else
                error('Not possible axisRatio with this volume')
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.ratio = cParams.superEllipseRatio;
            obj.volume = cParams.volumeMicro;
        end
        
        function isValid = isMicroStructureValid(obj,phi,volume)
            rMax = 0.99;
            phiMin = atan((volume)/(rMax^2));
            phiMax = atan((rMax^2)/(volume));
            isValid = phi <= phiMax && phi >= phiMin;
        end
        
    end
    
end