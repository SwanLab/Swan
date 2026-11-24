classdef DensityAspectSurrogate < handle
    properties (Access = private)
        model
        mU
    end
    methods (Access = public)
        function obj = DensityAspectSurrogate(model,mU)
            obj.model = model;
            obj.mU    = mU;
        end
        function Y = computeOutputValues(obj,X)
            Y = obj.model.computeOutputValues(obj.mapInput(X));
        end
        function [Y,dY] = computeOutputValuesAndGradient(obj,X)
            [Y,dYm] = obj.model.computeOutputValuesAndGradient(obj.mapInput(X));
            dY      = pagemtimes(obj.inputJacobian(X),'transpose',dYm,'none');
        end
        function dY = computeGradient(obj,X)
            [~,dY] = obj.computeOutputValuesAndGradient(X);
        end
        function E = computeError(obj,X,Y)
            E = norm(obj.computeOutputValues(X)-Y,'fro')/max(norm(Y,'fro'),eps);
        end
        function [m1,m2] = computeSides(obj,rho,t)
            [xi,~,~] = obj.computeXi(rho,t);
            m1 = sqrt((1-rho).*tan(xi));
            m2 = sqrt((1-rho)./tan(xi));
        end
    end
    methods (Access = private)
        function Xm = mapInput(obj,X)
            [m1,m2] = obj.computeSides(X(:,1),X(:,2));
            Xm      = [m1,m2,X(:,3:end)];
        end
        function J = inputJacobian(obj,X)
            rho     = X(:,1);
            t       = X(:,2);
            hole    = 1 - rho;
            [xi,dXi,range] = obj.computeXi(rho,t);
            tanXi   = tan(xi);
            slope   = hole.*(1 + tanXi.^2);
            [m1,m2] = obj.computeSides(rho,t);
            J       = repmat(eye(size(X,2)),1,1,size(X,1));
            J(1,1,:) = (-tanXi + slope.*dXi)./(2*m1);
            J(1,2,:) = slope.*range./(2*m1);
            J(2,1,:) = (-1./tanXi - slope.*dXi./tanXi.^2)./(2*m2);
            J(2,2,:) = -slope.*range./(2*m2.*tanXi.^2);
        end
        function [xi,dXi,range] = computeXi(obj,rho,t)
            hole  = 1 - rho;
            k     = obj.mU^2;
            xiMin = atan(hole/k);
            xiMax = atan(k./hole);
            range = xiMax - xiMin;
            xi    = xiMin + t.*range;
            dXi   = (2*t-1).*k./(k^2 + hole.^2);
        end
    end
end