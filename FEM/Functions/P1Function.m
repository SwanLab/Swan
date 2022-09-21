classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        type
        connec
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end

        function fC = computeValueInCenterElement(obj)
            % Goal: delete this function
            % Yields a different result from projecting to P0 (possibly due
            % to quadrature order?)
            q = Quadrature.set(obj.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fCenter = obj.evaluate(xV);
            fC = squeeze(fCenter);
        end

        function dF = computeDiscontinuousField(obj)
            % Goal: use this function
            s.fValues = obj.fValues;
            s.connec = obj.connec;
            s.type   = obj.type;
            dF = P1DiscontinuousFunction(s);
        end

        function plot(obj, m) % 2D domains only
            % Goal: use the P1DiscontinuousFunction instead
            dim = 1;
            figure()
            trisurf(m.connec, m.coord(:,1), m.coord(:,2), obj.fValues(:,dim))
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec = cParams.connec;
            obj.type   = cParams.type;
            obj.fValues = cParams.fValues;
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function fxV = evaluate(obj,xV) % Previously interpolateFunction
            % Note: function moved to P1DiscontinuousFunction
            % What shall we do here?
            fDfun = obj.computeDiscontinuousField();
            fVals = fDfun.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(fVals,1);
            nElem  = size(fVals,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = fVals(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end

    end

end