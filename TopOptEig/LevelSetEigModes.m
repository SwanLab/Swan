classdef LevelSetEigModes < DesignVariable
    
    properties (Access = private)
        creatorSettings
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = LevelSetEigModes(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.creatorSettings      = cParams.creatorSettings;
            obj.creatorSettings.type = cParams.initialCase;
            obj.createValue();
            obj.createUnfittedMesh();
        end
        
        function update(obj,value)
            obj.value = value;
            if ~isempty(obj.isFixed)
                obj.value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            obj.updateUnfittedMesh();
        end
        
        function m = getUnfittedMesh(obj)
            m = obj.unfittedMesh;
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value;
        end
        
        function vf = computeVolumeFraction(obj)
            q = Quadrature.set(obj.unfittedMesh.backgroundMesh.type);
            q.computeQuadrature('CONSTANT');
            bM = obj.unfittedMesh.backgroundMesh;
            dv = obj.unfittedMesh.computeDvolume(q);
            dVT = bM.computeDvolume(q)';
            vf = dv./dVT;
        end

        function gamma = getFirstEigenMode(obj)
            x = obj.value;
            gamma = x(end);
        end

    end
    
    methods (Access = private)
        
        function createValue(obj)
            s         = obj.creatorSettings;
            s.ndim    = obj.mesh.ndim;
            s.coord   = obj.mesh.coord;
            lsCreator = LevelSetCreator.create(s);
            lsValue   = lsCreator.getValue();
            obj.value = lsValue;
            obj.value(end+1,1) = 1;
        end
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh.innerMeshOLD;
            s.boundaryMesh   = obj.mesh.boxFaceMeshes;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.compute(obj.value);
        end
        
    end
    
end