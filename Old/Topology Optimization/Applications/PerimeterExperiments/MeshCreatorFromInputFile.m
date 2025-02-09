classdef MeshCreatorFromInputFile < handle

    properties (Access = public)
        backgroundMesh
        boundaryMesh
    end

    properties (Access = private)
        inputFile
    end

    methods (Access = public)

        function obj = MeshCreatorFromInputFile(cParams)
            obj.init(cParams);
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.inputFile = cParams.inputFile;
        end

        function createBackgroundMesh(obj)
            coord = [];
            connec = [];
            eval(obj.inputFile);
            s.coord  = coord(:,2:3);
            s.connec = connec(:,2:end);
            obj.backgroundMesh = Mesh.create(s);
        end

        function createBoundaryMesh(obj)
            eval(obj.inputFile);
            if exist('External_border_nodes','var')
                s.borderNodes    = External_border_nodes;
                s.borderElements = External_border_elements;
                s.backgroundMesh = obj.backgroundMesh;
                b = BoundaryMeshCreatorFromData(s);
                obj.boundaryMesh = b.create();
            else
                s.backgroundMesh = obj.backgroundMesh;
                s.dimension = 1:obj.backgroundMesh.ndim;
                bC = BoundaryMeshCreatorFromRectangularBox(s);
                obj.boundaryMesh = bC.create();
            end
        end



    end


end