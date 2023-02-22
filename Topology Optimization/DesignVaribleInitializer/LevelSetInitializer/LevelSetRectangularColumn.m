classdef LevelSetRectangularColumn < LevelSetCreator
    
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
            a = obj.value(1:nElem);
            b = obj.value(nElem+1:2*nElem);
            A = obj.computeVector(nElem,a);
            B = obj.computeVector(nElem,b);
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.computeRectangleCentre(x,y);
            ls1 = (x-obj.x0)-B/2;
            ls2 = (y-obj.y0)-A/2;
            ls3 = -(x-obj.x0)-B/2;
            ls4 = -(y-obj.y0)-A/2;
            ls = max([ls1,ls2,ls3,ls4],[],2);
            obj.levelSet = ls;
        end

        function computeRectangleCentre(obj,x,y)
            obj.x0 = 0.5*(max(x)+min(x));
            obj.y0 = 0.5*(max(y)+min(y));
        end

        function V = computeVector(obj,nElem,val)
            nnodes = size(obj.nodeCoord,1);
            V = zeros(nnodes,1);
            k = nnodes/nElem;
            c1 = 1;
            c2 = k;
            for i = 1:nElem
                V(c1:c2) = val(i);
                c1 = c1+k;
                c2 = c2+k;
            end
        end
    end
    
end