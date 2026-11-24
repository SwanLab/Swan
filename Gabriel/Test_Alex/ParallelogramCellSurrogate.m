classdef ParallelogramCellSurrogate < handle
    properties (Access = private)
        model
    end
    methods (Access = public)
        function obj = ParallelogramCellSurrogate(model)
            obj.model = model;
        end
        function Y = computeOutputValues(obj,X)
            Y = obj.model.computeOutputValues(obj.mapInput(X));
        end
        function [Y,dY] = computeOutputValuesAndGradient(obj,X)
            [Y,dYf] = obj.model.computeOutputValuesAndGradient(obj.mapInput(X));
            dY      = pagemtimes(obj.inputJacobian(X),'transpose',dYf,'none');
        end
        function dY = computeGradient(obj,X)
            [~,dY] = obj.computeOutputValuesAndGradient(X);
        end
        function E = computeError(obj,X,Y)
            E = norm(obj.computeOutputValues(X)-Y,'fro')/max(norm(Y,'fro'),eps);
        end
    end
    methods (Static, Access = public)
        function Xf = mapInput(X)
            a     = X(:,3);
            norm2 = (1-a).^2 + a.^2;
            Xf    = [X(:,1:2),(1-2*a)./norm2,2*a.*(1-a)./norm2,X(:,4),log(X(:,5))];
        end
    end
    methods (Static, Access = private)
        function J = inputJacobian(X)
            a     = X(:,3);
            norm2 = (1-a).^2 + a.^2;
            J     = zeros(6,5,size(X,1));
            J(1,1,:) = 1;
            J(2,2,:) = 1;
            J(3,3,:) = -4*a.*(1-a)./norm2.^2;
            J(4,3,:) = 2*(1-2*a)./norm2.^2;
            J(5,4,:) = 1;
            J(6,5,:) = 1./X(:,5);
        end
    end
end