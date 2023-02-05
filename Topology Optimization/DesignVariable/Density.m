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
            v{1} = obj.value;
        end
        
        function [fun, funNames] = getFunsToPlot(obj)
            aa.mesh = obj.mesh.meshes{1};
            aa.fValues = obj.value;
            valFun = P1Function(aa);

            fun = {valFun};
            funNames = {'Density'};
        end
        
        function rho = computeVolumeFraction(obj)
            s.mesh   = obj.mesh;
            s.fValues = obj.value;
            f = P1Function(s);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            rho = f.evaluate(xV);
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
                    obj.value  = 1 - heaviside(phi);
                case 'Given'
                    obj.value = s.rho0.*ones(size(obj.mesh.coord,1),1);
            end
        end
        
    end
    
end

