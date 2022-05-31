classdef DensityEigModes < DesignVariable
    
    properties (Access = private)
        creatorSettings
        initCase
    end
    
    methods (Access = public)
        
        function obj = DensityEigModes(cParams)
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

        function rho = getDensity(obj)
            x = obj.value;
            rho = x(1:end-1);
        end
        
        function gamma = getFirstEigenMode(obj)
            x = obj.value;
            gamma = x(end);
        end
        
    end
    
    methods (Access = protected)
        
       function init(obj,cParams)
            obj.type    = cParams.type;
            obj.mesh    = cParams.mesh;
            if isfield(cParams,'isFixed')            
              obj.isFixed = cParams.isFixed;
            end
            obj.initValue(cParams);
            if isprop(cParams,'scalarProductSettings')
                obj.createScalarProduct(cParams);
            end
       end

    end

    methods (Access = private)

        function initValue(obj,cParams)
            if isfield(cParams,'value')
                if isempty(cParams.value)
                    obj.value = ones(size(obj.mesh.coord,1)+1,1);
                else
                    obj.value = cParams.value;
                end
            else
                obj.value = ones(size(obj.mesh.coord,1)+1,1);
            end
        end

         function createScalarProduct(obj,cParams)
%             s = cParams.scalarProductSettings;
%             s.nVariables = obj.nVariables;
%             s.femSettings.mesh = obj.mesh;
%             obj.scalarProduct = ScalarProduct(s);
         end

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
                    obj.value = s.rho0.*ones(size(obj.mesh.coord,1)+1,1); % final number is the cost
            end
        end
        
    end
    
end