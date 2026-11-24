classdef PeriodicAngleSurrogate < handle
    properties (Access = private)
        model
        harmonic
    end
    methods (Access = public)
        function obj = PeriodicAngleSurrogate(model,harmonic)
            obj.model    = model;
            obj.harmonic = harmonic;
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
    methods (Access = private)
        function Xf = mapInput(obj,X)
            theta = X(:,end);
            Xf    = [X(:,1:end-1),cos(obj.harmonic*theta),sin(obj.harmonic*theta)];
        end
        function J = inputJacobian(obj,X)
            theta   = X(:,end);
            nInputs = size(X,2);
            J       = zeros(nInputs+1,nInputs,size(X,1));
            for i = 1:nInputs-1
                J(i,i,:) = 1;
            end
            J(nInputs,nInputs,:)   = -obj.harmonic*sin(obj.harmonic*theta);
            J(nInputs+1,nInputs,:) =  obj.harmonic*cos(obj.harmonic*theta);
        end
    end
end