classdef CharacteristicFunction < L2Function

    properties (Access = private)
        coorP1
    end

    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        mesh
        fxy
    end

    methods (Access = public)

        function obj = CharacteristicFunction(cParams)
            obj.init(cParams);
            obj.createP1CoorFunction();
        end

        function fxV = evaluate(obj,xV)
            xy    = obj.coorP1.evaluate(xV);
            nGaus = size(xy,2);
            nElem = size(xy,3);
            fxV   = zeros(1,nGaus,nElem);
            for iGaus = 1:nGaus
                x = xy(1,iGaus,:);
                y = xy(2,iGaus,:);
                f = obj.fxy(x,y);
                fxV(1,iGaus,f>0) = 0;
                fxV(1,iGaus,f<=0) = 1;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf = 1;
            obj.mesh  = cParams.mesh;
            obj.fxy   = cParams.fxy;
        end

        function createP1CoorFunction(obj)
            s.connec   = obj.mesh.connec;
            s.type     = obj.mesh.type;
            s.fValues  = obj.mesh.coord;
            s.ndimf    = size(obj.mesh.coord,2);
            obj.coorP1 = P1Function(s);
        end

    end
end