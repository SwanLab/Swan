classdef LevelSetWithSeveralHoles < LevelSetCreator
    
    properties (Access = private)
        hasToShowHoleInBCWarning
        bc
        nHoles
        rHoles
        phaseHoles
    end
    
    methods (Access = public)
        
        function obj = LevelSetWithSeveralHoles(cParams)
            obj.load_holes_settings(cParams);
            obj.loadWarningOption(cParams);
            obj.loadBoundaryConditions(cParams);
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeLevelSetValue()
            obj.showPossibleHoleinBcWarning()
        end
        
    end
    
    methods (Access = private)
        
        function loadBoundaryConditions(obj,cParams)
            s = cParams.geomParams;
            bCond = [];
            if ~isempty(s.dirichlet) && ~isempty(s.pointload)
                bCond = unique([s.dirichlet(:,1); s.pointload(:,1)]);
            end
            obj.bc = bCond;
        end
        
        function loadWarningOption(obj,cParams)
            s = cParams.geomParams;
            obj.hasToShowHoleInBCWarning = false;            
            if isfield(s,'warningHoleBC')
                if ~isempty(s.warningHoleBC)
                    obj.hasToShowHoleInBCWarning = s.warningHoleBC;
                end
            end
        end
        
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
        
        function showPossibleHoleinBcWarning(obj)
            if any(obj.levelSet(obj.bc)>0) && obj.hasToShowHoleInBCWarning
                warning('At least one BC is set on a hole')
            end
        end
        
        function load_holes_settings(obj,cParams)
            s = cParams.geomParams;
            obj.nHoles = s.nHoles;
            obj.rHoles = s.rHoles;
            obj.phaseHoles = s.phaseHoles;
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

