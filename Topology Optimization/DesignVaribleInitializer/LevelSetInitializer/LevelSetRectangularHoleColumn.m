classdef LevelSetRectangularHoleColumn < LevelSetCreator
    
    properties (Access = private)
        x0
        y0
    end

    properties (Access = private)
        radius
        value
    end
    
    methods (Access = public)
        
        function obj = LevelSetRectangularColumn(cParams)
            obj.value = cParams.desVarValue;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            nElem = size(obj.value,1)/2;
            h = obj.value(1:nElem);
            b = obj.value(nElem+1:2*nElem);
            eh = obj.value(2*nElem+1:3*nElem);
            eb = obj.value(3*nElem+1:4*nElem);
            outR = createOuterRectangleLevelSet(obj,h,b,nElem);
            inR  = -createInnerRectangleLevelSet(obj,eh,eb,nElem);
            ls = max([outR,inR],[],2);
            obj.levelSet = ls;
        end

        function levelSet = createOuterRectangleLevelSet(obj,h,b,nElem)
            s.type        = 'rectangularColumn';
            s.desVarValue = [h;b];
            s.coord       = obj.backgroundMesh.coord(1:2*nElem,1);
            s.ndim        = obj.backgroundMesh.ndim;
            lsCreator     = LevelSetCreator.create(s);
            levelSet      = lsCreator.getValue();            
        end
        
        function levelSet = createInnerRectangleLevelSet(obj,eh,eb,nElem)
            s.type        = 'rectangularColumn';
            s.desVarValue = [eh;eb];
            s.coord       = obj.backgroundMesh.coord(2*nElem+1:4*nElem,1);
            s.ndim        = obj.backgroundMesh.ndim;
            lsCreator     = LevelSetCreator.create(s);
            levelSet      = lsCreator.getValue();            
        end
    end
    
end