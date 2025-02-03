classdef FinalMeshTestComputer < handle

    properties (Access = public)
        comparingVerticesMessage
        comparingFaceMessage
    end

    properties (Access = private)
        tolerance
    end

    properties (Access = private)
        finalMeshOriginal
        finalMeshTest
    end
    
    methods (Access = public)

        function obj = FinalMeshTestComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.compareVertices();
            obj.compareFaces();
            result = [obj.comparingVerticesMessage, newline, obj.comparingFaceMessage];
            disp(result);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.saveInput(cParams);
            obj.loadOriginalfinalMesh();
            obj.tolerance = 1e-6;
        end

        function saveInput(obj,cParams)
            obj.finalMeshTest    = cParams.mesh;
        end

        function loadOriginalfinalMesh(obj)
            load("datas.mat","mesh");
            obj.finalMeshOriginal = mesh;
        end

        function compareVertices(obj)
            try
                assert(all(abs(obj.finalMeshOriginal.coord - obj.finalMeshTest.coord) < obj.tolerance, 'all'), ...
                'The final mesh does not have the expected value for vertex positions.');
                 obj.comparingVerticesMessage = 'The final mesh has the expected value for vertex positions.';
            catch ME
                obj.comparingVerticesMessage = ME.message;
            end
        end

        function compareFaces(obj)
            try
                assert(all(abs(obj.finalMeshOriginal.connec - obj.finalMeshTest.connec) < obj.tolerance, 'all'), ...
                    'The final mesh does not have the expected value for face connectivity.');
                obj.comparingFaceMessage = 'The final mesh has the expected value for face connectivity.';
            catch ME
                obj.comparingFaceMessage = ME.message;
            end
        end

    end
    
end