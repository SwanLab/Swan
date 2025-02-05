classdef ReferenceMeshTestComputer < handle

    properties (Access = public)
        comparingVerticesMessage
        comparingFaceMessage
    end

    properties (Access = private)
        tolerance
    end

    properties (Access = private)
        referenceMeshOriginal
        referenceMeshTest
    end
    
    methods (Access = public)

        function obj = ReferenceMeshTestComputer(cParams)
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
            obj.loadOriginalReferenceMesh();
            obj.tolerance = 1e-6;
        end

        function saveInput(obj,cParams)
            obj.referenceMeshTest    = cParams.refMesh;
        end

        function loadOriginalReferenceMesh(obj)
            load("datas.mat","m");
            obj.referenceMeshOriginal = m;
        end

        function compareVertices(obj)
            try
                assert(all(abs(obj.referenceMeshOriginal.coord - obj.referenceMeshTest.coord) < obj.tolerance, 'all'), ...
                'The reference mesh does not have the expected value for vertex positions.');
                 obj.comparingVerticesMessage = 'The reference mesh has the expected value for vertex positions.';
            catch ME
                obj.comparingVerticesMessage = ME.message;
            end
        end

        function compareFaces(obj)
            try
                assert(all(abs(obj.referenceMeshOriginal.connec - obj.referenceMeshTest.connec) < obj.tolerance, 'all'), ...
                    'The reference mesh does not have the expected value for face connectivity.');
                obj.comparingFaceMessage = 'The reference mesh has the expected value for face connectivity.';
            catch ME
                obj.comparingFaceMessage = ME.message;
            end
        end

    end
    
end