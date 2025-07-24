classdef LevelSetWithSeveralHoles < LevelSetCreator
    
    properties (Access = private)
        nHoles
        rHoles
        phaseHoles
    end
    
    methods (Access = public)
        
        function obj = LevelSetWithSeveralHoles(cParams)
            obj.init(cParams);
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeLevelSetValue();
        end
        
    end
    
    methods (Access = private)
        
        function computeLevelSetValue(obj)
            ls = ones(obj.lsSize);
            for idim = 1:obj.ndim
                coordV = obj.nodeCoord(:,idim);
                cosDim = obj.computeDirectionalCosinus(coordV,idim);
                ls = ls.*cosDim;
            end
            ls = ls + obj.rHoles-1;
            obj.levelSet = ls;
        end
        
        function init(obj,cParams)
            obj.nHoles     = cParams.nHoles;
            obj.rHoles     = cParams.rHoles;
            obj.phaseHoles = cParams.phaseHoles;
        end
        
        function cosDir  = computeDirectionalCosinus(obj,coord,dir)
            pos = coord;
            l = max(pos) - min(pos);
            n = obj.nHoles(dir);
            fase = obj.phaseHoles(dir);
            cosDir = cos((n + 1)*pos*pi/l + fase);
        end
        
    end
    
end

