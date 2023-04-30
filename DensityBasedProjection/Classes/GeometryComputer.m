classdef GeometryComputer < handle
    properties (Access = public)
        degress
        connectivityDOFS
        elementNumberX
        elementNumberY
    end
    properties (Access = private)
        nodeToDegrees
        data
        elementType
        connectivityMesh
        boundaryConditions
        coord
        
    end

    methods (Access = public)
        function obj = GeometryComputer(cParams)
            obj.inputData(cParams);
            obj.createNodeToDegreesMat();
        end
        function compute(obj)
            obj.countElements();
            obj.connectDegrees();
            obj.computeDregreess();
        end
    end
    methods (Access = private)
        function obj = inputData(obj,cParams)
            obj.connectivityMesh = cParams.mesh.connec;
            obj.coord = cParams.mesh.coord;
            obj.boundaryConditions = cParams.bc;
            obj.elementType = cParams.mesh.type;
        end
        function createNodeToDegreesMat(obj)
            linearMat = (1:1:length(obj.data.gidcoord)*2);
            obj.nodeToDegrees = reshape(linearMat,2,[]);
        end
        function countElements(obj)
            obj.elementNumberX = length(unique(obj.coord(:,1)))-1;
            obj.elementNumberY = length(unique(obj.coord(:,2)))-1;

        end 
        function computeDregreess(obj)

            obj.degress.fixed = reshape(obj.nodeToDegrees(:,obj.data.nodesolid),1,[]);
            obj.degress.all   =  reshape(obj.nodeToDegrees,1,[]);
            obj.degress.free  = setdiff(obj.degress.all,obj.degress.fixed);

        end
        function connectDegrees(obj)
            switch obj.elementType
                case 'QUAD'
                    obj.connectDegreesSquare();
                case 'TRIANGLE'
                    obj.connectDegreesTriangular();
                otherwise
                    disp('Element mesh not implemented');
            end
        end
        function connectDegreesSquare(obj)
            for e = 1:length(obj.data.gidlnods)
                obj.connectivityDOFS(e,:) = reshape(obj.nodeToDegrees(:,obj.connectivityMesh(e,:)),1,[]);
            end
        end
        function connectDegreesTriangular(obj)
            for e = 1:length(obj.data.gidlnods)
                obj.connectivityDOFS(e,:) = reshape(obj.nodeToDegrees(:,obj.connectivityMesh(e,:)),1,[]);
            end
        end
    end
end
