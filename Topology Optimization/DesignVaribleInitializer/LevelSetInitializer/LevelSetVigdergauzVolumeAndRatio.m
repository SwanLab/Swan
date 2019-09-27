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
                x = obj.nodeCoord(:,1);                
                y = obj.nodeCoord(:,2);                
                obj.levelSet = LevelSetLipung(x,y,volum,phi);
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