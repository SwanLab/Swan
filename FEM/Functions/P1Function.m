classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        type
        connec

        fNodes %just to check computeDiscontField
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createDiscontinuousP1();
        end

        function fC = computeValueInCenterElement(obj)
            % Goal: delete this function
            q = Quadrature.set(obj.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fCenter = obj.evaluate(xV);
            fC = squeeze(fCenter);
        end

        function fxV = evaluate(obj,xV) % Previously interpolateFunction
            % Note: the same function is copied at P1DiscontinuousFunction
            % What shall we do here?
            func = obj.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end

        function dF = computeDiscontinuousField(obj)
            % Goal: use this function
            s.fNodes = obj.fNodes;
            s.connec = obj.connec;
            dF = P1DiscontinuousFunction(s);
        end

        function plot(obj, m) % 2D domains only
            % Goal: use the P1DiscontinuousFunction instead
            dim = 1;
            figure()
            trisurf(m.connec, m.coord(:,1), m.coord(:,2), obj.fNodes(:,dim))
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec = cParams.connec;
            obj.type   = cParams.type;
            obj.fValues = cParams.fNodes;
            obj.fNodes = cParams.fNodes;
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function createDiscontinuousP1(obj)
            % Goal: merge with computeDiscontinuousField
            % Goal: the fValues of this class should be the current fNodes
            s.fNodes = obj.fNodes;
            s.connec = obj.connec;
            p1d = P1DiscontinuousFunction(s);
            obj.fValues = p1d.fValues;
        end

    end

end