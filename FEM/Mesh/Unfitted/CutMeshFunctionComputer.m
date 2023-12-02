classdef CutMeshFunctionComputer < handle

    properties (Access = private)
        backgroundMesh
        cutMesh
    end

    properties (Access = private)
        coordGlobal
    end

    methods (Access = public)
        function obj = CutMeshFunctionComputer(cParams)
            obj.init(cParams);
            obj.computeGlobalCoordinatesMatrix();
            % ...
        end

        function newFeFun = compute(obj,backgroundFeFun)
            oldfValues = backgroundFeFun.fValues;
            % ...
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutMesh        = cParams.cutMesh;
        end

        function computeGlobalCoordinatesMatrix(obj)
            coord           = obj.backgroundMesh.coord;
            icCoord         = obj.cutMesh.coord;
            ind             = not(ismember(icCoord,coord,'rows'));
            obj.coordGlobal = [coord;icCoord(ind,:)];
        end
    end
end