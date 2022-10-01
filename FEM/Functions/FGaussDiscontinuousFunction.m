classdef FGaussDiscontinuousFunction < handle
    
    properties (Constant, Access = public)
        fType = 'GAUSSPOINTS'
    end

    properties (Access = public)
        ndimf
        fValues
        quadrature
    end

    properties (Access = private)
        type % not needed apparently
        connec
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = FGaussDiscontinuousFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xV)
            assert(isequal(xV, obj.quadrature.posgp), 'Gauss points do not match')
            fxV = obj.fValues;
        end

        function plot(obj, m)
            % Plot using a P1 Disc function
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues    = cParams.fValues;
            obj.connec     = cParams.connec;
            obj.type       = cParams.type;
            obj.quadrature = cParams.quadrature;
            obj.ndimf      = size(cParams.fValues,1);
        end
        
    end
    
end