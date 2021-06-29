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
        
        function rho = computeVolumeFraction(obj)
            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = obj.value;
            f = FeFunction(s);
            rho = f.computeValueInCenterElement();            
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

