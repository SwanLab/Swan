classdef TestPlotting < handle
    
    properties (Access = private)
        testName;
        meshType;
        meshIncludeBoxContour;
    end

    properties (Access = private)
        varAdim
        mesh
        settings
        computation
        levelSet
        unfittedMesh
        oldMeshUnfitted
    end
    
    methods (Access = public)

        
        function obj = TestPlotting(cParams)
            obj.init(cParams);
            obj.createLevelSet();
            obj.createMesh();
            obj.plot();
        end
        
        function plot(obj)
            figure();
            if isequal(obj.meshType,'BOUNDARY')
                obj.unfittedMesh.plotBoundary();
                obj.unfittedMesh.plotNormals();
            elseif isequal(obj.meshType,'INTERIOR')
                obj.unfittedMesh.plot();
            end
            view(obj.getViewAngle());
        end

        function passed = computePassed(obj)
            d = load(obj.testName);
            unfittedMesh = obj.unfittedMesh;
            loaded = d.unfittedMesh;
            passed = isequaln(unfittedMesh,loaded);
            if ~passed
                save(obj.testName, 'unfittedMesh', '-append')
            end
        end

        function overwriteResults(obj)
            save(obj.testName, 'unfittedMesh', '-append')
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.testName        = cParams.testName;
            obj.meshType        = cParams.meshType;
            obj.meshIncludeBoxContour = cParams.meshIncludeBoxContour;
        end

        function createLevelSet(obj)
            run(obj.testName)
            
            fun = GeometricalFunction(geomFunSettings);
            lsfun = fun.computeLevelSetFunction(mMesh);
            obj.mesh = mMesh;
            obj.levelSet = lsfun.fValues;
        end

        function createMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uM = UnfittedMesh(s);
            uM.compute(obj.levelSet);
            obj.unfittedMesh = uM;
        end

        function angle = getViewAngle(obj)
            if isprop(obj,'viewAngle')
                angle = obj.viewAngle;
            else
                angle = [0 0 1];
            end
        end
        
    end

end