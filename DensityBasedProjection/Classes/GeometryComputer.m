classdef GeometryComputer < handle
    properties (Access = public)
        conectivityMatrixMat
        degress
    end
    properties (Access = private)
        elementNumberX
        elementNumberY       
    end

    methods (Access = public)
        function obj = GeometryComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.connectNodes();
            obj.computeDregreess();
        end
    end
    methods (Access = private)
        function obj = inputData(obj,cParams)
            obj.elementNumberX =  cParams.elementNumberX;
            obj.elementNumberY =  cParams.elementNumberY;
        end
        function connectNodes(obj)
            nodeNumeration = reshape(1:(1+obj.elementNumberX)*(1+obj.elementNumberY),1+obj.elementNumberY,1+obj.elementNumberX);
            conectivityMatrixVec = reshape(2*nodeNumeration(1:end-1,1:end-1)+1,obj.elementNumberX*obj.elementNumberY,1);
            obj.conectivityMatrixMat = repmat(conectivityMatrixVec,1,8)+repmat([0 1 2*obj.elementNumberY+[2 3 0 1] -2 -1],obj.elementNumberX*obj.elementNumberY,1);
        end 
        function computeDregreess(obj)
            obj.degress.fixed = [1:2:2*(obj.elementNumberY+1) 2*(obj.elementNumberX+1)*(obj.elementNumberY+1) 2*(obj.elementNumberX+1)*(obj.elementNumberY+1)-1];
            obj.degress.all   = 1:2*(obj.elementNumberX+1)*(obj.elementNumberY+1);
            obj.degress.free  = setdiff(obj.degress.all,obj.degress.fixed);
            
        end 
    end
end