classdef GeometryComputer < handle
    properties (Access = public)
        degress
        connectivityDOFS
        elementNumberX
        elementNumberY
        elementType
    end
    properties (Access = private)
        nodeToDegrees
        data
        
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
            linearMat = (1:1:length(obj.coord)*2);
            obj.nodeToDegrees = reshape(linearMat,2,[]);
        end
        function countElements(obj)
            obj.elementNumberX = length(unique(obj.coord(:,1)))-1;
            obj.elementNumberY = length(unique(obj.coord(:,2)))-1;

        end 
        function computeDregreess(obj)
            
            obj.degress.all   =  reshape(obj.nodeToDegrees,1,[]);
            for e = 1:length(obj.boundaryConditions.dirichlet)
            obj.degress.fixed(e) = obj.nodeToDegrees(obj.boundaryConditions.dirichlet(e,2),obj.boundaryConditions.dirichlet(e,1));
            end 
            obj.degress.free  = setdiff(obj.degress.all,obj.degress.fixed);

            for e = 1:length(obj.boundaryConditions.pointload)
            obj.degress.forceDOFs(e) = obj.nodeToDegrees(obj.boundaryConditions.pointload(e,2),obj.boundaryConditions.pointload(e,1));
            end 


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
            for e = 1:length(obj.connectivityMesh)
                obj.connectivityDOFS(e,:) = reshape(obj.nodeToDegrees(:,obj.connectivityMesh(e,:)),1,[]);
            end
        end
        function connectDegreesTriangular(obj)
            for e = 1:length(obj.connectivityMesh)
                obj.connectivityDOFS(e,:) = reshape(obj.nodeToDegrees(:,obj.connectivityMesh(e,:)),1,[]);
            end
        end
    end
end
