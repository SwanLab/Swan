classdef TestVigdergauzMicroStructure < handle
    
    properties (Access = private)
        designVariable
        unfittedMesh
        volume
        numericalVolume
        computation
        settings
        testName;
    end
    
    methods (Access = public)
        
        function obj = TestVigdergauzMicroStructure(cParams)
            obj.init(cParams);
            obj.createUnfittedMesh();
            obj.computeFractionVolume();
        end

        function error = computeError(obj)
            v = obj.volume;
            nv = obj.numericalVolume;
            error = abs(v - nv)/v;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.volume   = cParams.volume;
            obj.testName = cParams.testName;
            s.testName = obj.testName;
            computer = TopOptComputer(s);
            computer.compute();
            obj.computation = computer.computation;
            obj.settings = computer.settings;
            obj.designVariable = obj.computation.designVariable;
        end

        function createUnfittedMesh(obj)
            meshBackground = obj.computation.designVariable.mesh;
            s.backgroundMesh  = meshBackground.innerMeshOLD;
            s.boundaryMesh    = meshBackground.boxFaceMeshes;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.unfittedMesh.compute(obj.designVariable.value); % peta aqui
        end

        function computeFractionVolume(obj)
            v = obj.unfittedMesh.computeMass();
            obj.numericalVolume = v;
        end

    end

end