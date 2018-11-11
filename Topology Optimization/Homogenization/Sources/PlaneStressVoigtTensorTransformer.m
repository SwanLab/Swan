classdef PlaneStressVoigtTensorTransformer < handle
    
    properties (Access = public)
        tensorVoigtInPlaneStress
    end
    
    properties (Access = private)
       C 
       InPlane
       OutPlane
       nstre
       detOutPlane
       deter
    end
    
    
    methods (Access = public)
        function obj = PlaneStressVoigtTensorTransformer(C)
            obj.init(C)
            obj.computeOutPlaneDeterminant()
            obj.computeDeterminants()
            obj.computeComponents()
        end
    end
    
    methods (Access = private)
        
        function init(obj,C)
            obj.C = C;
            obj.InPlane  = [1 2 6];
            obj.OutPlane = [3 4 5];
            obj.nstre    = 3;
        end
        
        function computeOutPlaneDeterminant(obj)
            obj.detOutPlane = obj.compute3x3Determinant(obj.C(obj.OutPlane,obj.OutPlane));
        end
        
        function computeDeterminants(obj)
            obj.deter = sym(zeros(9,1));
              for i = 1:obj.nstre
                for j = 1:obj.nstre
                    InPlaneIndex  = obj.InPlane(i);
                    OutPlaneIndex = setdiff(obj.OutPlane,obj.OutPlane(j));
                    ColumnIndex   = [InPlaneIndex OutPlaneIndex];
                    RowIndex      = obj.OutPlane;
                    SubTensor     = obj.C(RowIndex,ColumnIndex);
                    StoredIndex   = obj.nstre*(j-1) + i;
                    det = obj.compute3x3Determinant(SubTensor);
                    obj.deter(StoredIndex) = det;
                end
              end
       end
        
        function component = computeComponent(obj,i,j)
            A   = sym(zeros(4,1));
            det = sym(zeros(4,1));
            for iDet = 1:4
                det(iDet) = obj.getDeterminant(iDet,j);
                A(iDet)   = obj.getTensorTerm(iDet,i,j);
            end
            component = obj.compute4x4Determinant(A,det);
        end
        
        function TensorTerm = getTensorTerm(obj,iDet,i,j)
               ColumnIndex = obj.computeColumnIndex(iDet,j);
               RowsIndex   = obj.InPlane(i);
               TensorTerm  = obj.C(RowsIndex,ColumnIndex);
        end
        
        function Index = computeColumnIndex(obj,iDet,j)
            if obj.isFirst(iDet)
                Index  = obj.InPlane(j);
            else
                Index = obj.OutPlane(iDet-1);
            end
        end
        
        function det = getDeterminant(obj,iDet,j)
            if obj.isFirst(iDet)
                det = obj.detOutPlane;
            else
                StoredIndex = obj.nstre*(iDet-2) + j;
                det = obj.deter(StoredIndex);
            end
        end
        
        function  computeComponents(obj)
            A = sym(zeros(obj.nstre,obj.nstre));
            for i = 1:obj.nstre
                for j=i:obj.nstre
                    A(i,j) = obj.computeComponent(i,j);
                    A(j,i) = A(i,j);
                end
            end
            A = A/(obj.detOutPlane);
            obj.tensorVoigtInPlaneStress = A;
        end
        
    end
    
    methods (Static, Access = private)
        function value = compute4x4Determinant(A,det)
            value = 0;
            for i = 1:length(det)
                value = value + (-1)^(i-1)*A(i)*det(i);
            end
        end
        
        function  det = compute3x3Determinant(A)
            det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)+ A(3,2)*A(2,1)*A(1,3)...
                 -A(1,3)*A(2,2)*A(3,1) - A(2,3)*A(3,2)*A(1,1)- A(2,1)*A(1,2)*A(3,3);
        end
        
        function isFirst = isFirst(index)
            isFirst = index == 1;
        end
    end
    
end

