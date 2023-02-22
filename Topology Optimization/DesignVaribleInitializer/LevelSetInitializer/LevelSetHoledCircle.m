classdef LevelSetHoledCircle < LevelSetCreator
    
     properties (Access = private)
        x0
        y0
    end

    properties (Access = private)
        radius
        value
    end
    
    methods (Access = public)
        
        function obj = LevelSetHoledCircle(cParams)
            obj.value = cParams.desVarValue;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            nElem = size(obj.value,1)/2;
            r1 = obj.value(1:nElem);
            r2 = r1 + obj.value(nElem+1:2*nElem);
            R1 = obj.computeRadiusVector(nElem,r1);
            R2 = obj.computeRadiusVector(nElem,r2);
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.computeCircumferenceCentre(x,y);
            ls1 = (x-obj.x0).^2 + (y-obj.y0).^2 - R1.^2;
            ls1 = -ls1;
            ls2 = (x-obj.x0).^2 + (y-obj.y0).^2 - R2.^2;
            ls = max([ls1,ls2],[],2);
            obj.levelSet = ls;
        end

        function computeCircumferenceCentre(obj,x,y)
            obj.x0 = 0.5*(max(x)+min(x));
            obj.y0 = 0.5*(max(y)+min(y));
        end

        function R = computeRadiusVector(obj,nElem,radius)
            nnodes = size(obj.nodeCoord,1);
            R = zeros(nnodes,1);
            k = nnodes/nElem;
            c1 = 1;
            c2 = k;
            for i = 1:nElem
                R(c1:c2) = radius(i);
                c1 = c1+k;
                c2 = c2+k;
            end
        end
    end
    
end