classdef Density < DesignVariable
    
    properties (Access = private)
        creatorSettings
        initCase
    end
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.initCase = cParams.initialCase;
            obj.creatorSettings  = cParams.creatorSettings;
            obj.createValue();
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
    
    methods (Access = private)
        
        function createValue(obj)
            s = obj.creatorSettings;
            switch s.type 
                case 'FromLevelSet'
                    s.ndim  = obj.mesh.ndim;
                    s.coord = obj.mesh.coord;
                    s.type  = obj.initCase;
                    lsCreator  = LevelSetCreator.create(s);
                    phi        = lsCreator.getValue();
                    value  = 1 - heaviside(phi);
                case 'Given'
                    value = s.rho0.*ones(size(obj.mesh.coord,1),1);
            end
            ss.mesh    = obj.mesh;
            ss.fValues = value;
            obj.fun    = P1Function(ss);
        end

    end
    
end

