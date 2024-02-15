classdef Density < DesignVariable
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.fun.fValues;
        end
        
        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.fun};
            funNames = {'Density'};
        end
        
        function rho = computeVolumeFraction(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            rho = obj.fun.evaluate(xV);
        end
        
    end
    
end

