classdef GeometryComputer < handle
    properties (Access = public)
        conectivityMatrixMat
        degress
    end
    properties (Access = private)
        elementNumberX
        elementNumberY    
        nodeToDegrees
        data
        elementType
    end

    methods (Access = public)
        function obj = GeometryComputer(cParams)
            obj.inputData(cParams);
            obj.createNodeToDegreesMat();
        end
        function compute(obj)
            obj.connectDegrees();
            obj.computeDregreess();
        end
    end
    methods (Access = private)
        function obj = inputData(obj,cParams)
            obj.data = cParams.data;
            obj.elementNumberX =  length(unique(obj.data.gidcoord(:,2)));
            obj.elementNumberY =  length(unique(obj.data.gidcoord(:,3)));
            obj.elementType = cell2mat(obj.data.Data_prb(1));
        end
        function createNodeToDegreesMat(obj)
            linearMat = (1:1:length(obj.data.gidcoord)*2);
            obj.nodeToDegrees = reshape(linearMat,2,[]);
        end 
%         function connectNodes(obj)
%             nodeNumeration = reshape(1:(1+obj.elementNumberX)*(1+obj.elementNumberY),1+obj.elementNumberY,1+obj.elementNumberX);
%             conectivityMatrixVec = reshape(2*nodeNumeration(1:end-1,1:end-1)+1,obj.elementNumberX*obj.elementNumberY,1);
%             obj.conectivityMatrixMat = repmat(conectivityMatrixVec,1,8)+repmat([0 1 2*obj.elementNumberY+[2 3 0 1] -2 -1],obj.elementNumberX*obj.elementNumberY,1);
%         end 

        function computeDregreess(obj)
            obj.degress.fixed = [1:2:2*(obj.elementNumberY+1) 2*(obj.elementNumberX+1)*(obj.elementNumberY+1) 2*(obj.elementNumberX+1)*(obj.elementNumberY+1)-1];
            obj.degress.all   = 1:2*(obj.elementNumberX+1)*(obj.elementNumberY+1);
            obj.degress.free  = setdiff(obj.degress.all,obj.degress.fixed);
        end
        function connectDegrees(obj)
            switch obj.elementType
                case 'SQUARE'
                    obj.connectDegreesSquare();
                case 'TRIANGLE'
                    obj.connectDegreesTriangular();
                otherwise
                    disp('Element mesh not implemented');
            end
        end 
        function connectDegreesSquare(obj)
            for e = 1:length(obj.data.gidlnods)
             obj.conectivityMatrixMat(e,:) = reshape(obj.nodeToDegrees(:,obj.data.gidlnods(e,2:end)),1,[]);
            end 
        end 
        function connectDegreesTriangular(obj)
            for e = 1:length(obj.data.gidlnods)
              obj.conectivityMatrixMat(e,:) = reshape(obj.nodeToDegrees(:,obj.data.gidlnods(e,2:end-1)),1,[]);
            end 
        end         
    end
end