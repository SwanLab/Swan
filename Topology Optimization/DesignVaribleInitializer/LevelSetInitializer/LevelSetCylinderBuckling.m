classdef LevelSetCylinderBuckling < LevelSetCreator
    
    properties (Access = private)
        x0
        y0
    end

    properties (Access = private)
        radius
        value
    end
    
    methods (Access = public)
        
        function obj = LevelSetCylinderBuckling(cParams)
            obj.value = cParams.desVarValue;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            nElem = size(obj.value,1);
            R = obj.computeRadiusVector(nElem);
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.computeCircumferenceCentre(x,y);
            ls = (x-obj.x0).^2 + (y-obj.y0).^2 - R.^2;
            obj.levelSet = ls;
        end

        function computeCircumferenceCentre(obj,x,y)
            obj.x0 = 0.5*(max(x)+min(x));
            obj.y0 = 0.5*(max(y)+min(y));
        end

        function R = computeRadiusVector(obj,nElem)
            nnodes = size(obj.nodeCoord,1);
            R = zeros(nnodes,1);
            k = nnodes/nElem;
            c1 = 1;
            c2 = k;
            for i = 1:nElem
                R(c1:c2) = obj.value(i);
                c1 = c1+k;
                c2 = c2+k;
            end
        end
    end
    
end