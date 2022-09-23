classdef FGaussDiscontinuousFunction < FeFunction
    
    properties (Access = private)
        connec
        type
    end
    
    properties (Access = private)
        xG
    end
    
    methods (Access = public)
        
        function obj = FGaussDiscontinuousFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xV)
            assert(xV==obj.quadrature.posgp, 'Gauss points do not match')
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
        end
        
    end
    
end