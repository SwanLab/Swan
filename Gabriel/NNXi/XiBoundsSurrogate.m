classdef XiBoundsSurrogate < handle

    properties (Access = private)
        model
        mU
    end

    methods (Access = public)

        function obj = XiBoundsSurrogate(model,mU)
            obj.model = model;
            obj.mU    = mU;
        end

        function Y = computeOutputValues(obj,X)
            Y = obj.model.computeOutputValues(obj.mapInput(X));
        end

        function [Y,dY] = computeOutputValuesAndGradient(obj,X)
            [Y,dYxi] = obj.model.computeOutputValuesAndGradient(obj.mapInput(X));
            dY = pagemtimes(obj.inputJacobian(X),'transpose',dYxi,'none');
        end

        function dY = computeGradient(obj,X)
            [~,dY] = obj.computeOutputValuesAndGradient(X);
        end

        function E = computeError(obj,X,Y)
            E = norm(obj.computeOutputValues(X)-Y,'fro')/max(norm(Y,'fro'),eps);
        end

        function xi = computeXi(obj,rho,t)
            [xiMin,xiMax] = obj.computeBounds(rho);
            xi = xiMin + t.*(xiMax-xiMin);
        end

    end

    methods (Access = private)

        function Xxi = mapInput(obj,X)
            Xxi        = X;
            Xxi(:,end) = obj.computeXi(X(:,1),X(:,end));
        end

        function J = inputJacobian(obj,X)
            rho = X(:,1);
            t   = X(:,end);
            nIn = size(X,2);
            [xiMin,xiMax]   = obj.computeBounds(rho);
            [dXiMin,dXiMax] = obj.computeBoundsDerivative(rho);
            J = repmat(eye(nIn),1,1,size(X,1));
            J(end,1,:)   = (1-t).*dXiMin + t.*dXiMax;
            J(end,end,:) = xiMax - xiMin;
        end

        function [xiMin,xiMax] = computeBounds(obj,rho)
            holeArea = 1 - rho;
            xiMin    = atan(holeArea/obj.mU^2);
            xiMax    = atan(obj.mU^2./holeArea);
        end

        function [dXiMin,dXiMax] = computeBoundsDerivative(obj,rho)
            holeArea = 1 - rho;
            m2       = obj.mU^2;
            dXiMin   = -m2./(m2^2 + holeArea.^2);
            dXiMax   =  m2./(holeArea.^2 + m2^2);
        end

    end

end