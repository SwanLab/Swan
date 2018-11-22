classdef LevelSetWithSeveralHoles < LevelSetCreator
    
    properties (Access = private)
        hasToShowHoleInBCWarning
        nHoles
        rHoles
        phaseHoles       
    end
    
    methods (Access = public)
        
        function obj = LevelSetWithSeveralHoles(input)
            obj.load_holes_settings(input);
            obj.loadWarningOption(input);
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet();
            obj.computeDesignVariable();
            obj.showPossibleHoleinBcWarning();
        end
        
    end
    
    methods (Access = private)
        
        function loadWarningOption(obj,input)
            if isempty(input.warningHoleBC)
                obj.hasToShowHoleInBCWarning = true;
            else
                obj.hasToShowHoleInBCWarning = input.warningHoleBC;
            end
        end
        
        function load_holes_settings(obj,input)
            obj.nHoles = input.nHoles;
            obj.rHoles = input.rHoles;
            obj.phaseHoles = input.phaseHoles;
        end
        

        function computeLevelSet(obj)
            ls = ones(obj.lsSize); 
            for idim = 1:obj.ndim
                coordV = obj.nodeCoord(:,idim);
                cosDim = obj.computeDirectionalCosinus(coordV,idim);
                ls = ls.*cosDim;
            end
            ls = ls + obj.rHoles-1;
            obj.levelSet = ls;
        end
        
        function cosDir  = computeDirectionalCosinus(obj,coord,dir)
            pos = coord;
            l = max(pos) - min(pos);
            n = obj.nHoles(dir);
            fase = obj.phaseHoles(dir);
            cosDir = cos((n + 1)*pos*pi/l + fase);
        end
              
        function showPossibleHoleinBcWarning(obj)
            bc = unique([obj.mesh.dirichlet(:,1); obj.mesh.pointload(:,1)]);
            if any(obj.x(bc)>0) && obj.hasToShowHoleInBCWarning
                warning('At least one BC is set on a hole')
            end
        end
    end
    
end

