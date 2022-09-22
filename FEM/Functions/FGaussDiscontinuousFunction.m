classdef FGaussDiscontinuousFunction < FeFunction
    
    properties (Access = private)
        connec
        type
        xG
    end
    
    methods (Access = public)
        
        function obj = FGaussDiscontinuousFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xV)
            obj.checkGaussPoints();
            nGaus  = size(obj.xG,2);
            func   = obj.fValues;
            for kGaus = 1:nGaus
                if xV == obj.xG(:,kGaus)
                    fxV = func(:,kGaus,:);
                end
            end
        end

%         function plot(obj, m)
            % I DON'T KNOW AT WHICH EXTENT IT MAKES SENSE TO PLOT A FGAUSSDISCFUN
%         end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
        end

        function checkGaussPoints(obj)
            q = Quadrature.set(obj.type);
            q.computeQuadrature('LINEAR');
            obj.xG = q.posgp;
        end
        
    end
    
end