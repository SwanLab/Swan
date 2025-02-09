classdef PST4VoigtFourthOrderTensorNumerically < PST4VoigtFourthOrderTensor
    
   properties (Access = private)
       C 
       InPlane
       OutPlane
       nstre
       detOutPlane
       deter
       isSymmbolic
    end
    
    
    methods (Access = public)
         
        function obj = PST4VoigtFourthOrderTensorNumerically(Tensor)
            obj.createPlaneStressTensor();
            obj.init(Tensor)
            obj.computeOutPlaneDeterminant()
            obj.computeDeterminants()
            obj.computeComponents()
        end
        
        function Index = getInPlaneIndex(obj)
            Index = obj.InPlane;
        end
    end
    
    
    methods (Access = private)
        
         function init(obj,tensor)
            PSIndex         = PlaneStressIndex();
            obj.InPlane     = PSIndex.getInPlaneIndex();
            obj.OutPlane    = PSIndex.getOutPlaneIndex();            
            tensSize        = obj.psTensor.getTensorSize();
            obj.isSymmbolic = isa(tensor,'sym');
            obj.nstre       = tensSize(1);
            obj.C = tensor;
        end
        
        function computeOutPlaneDeterminant(obj)
            obj.detOutPlane = obj.compute3x3Determinant(obj.C(obj.OutPlane,obj.OutPlane));
        end
        
        function computeDeterminants(obj)
            d = obj.nstre*obj.nstre;
            obj.deter = obj.createVector(d);
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
            d = 4;
            A   = obj.createVector(d);
            det = obj.createVector(d);
            for iDet = 1:d
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
            d = obj.nstre;
            A = obj.createSquareMatrix(d);
            for i = 1:d
                for j=i:d
                    A(i,j) = obj.computeComponent(i,j);
                    A(j,i) = A(i,j);
                end
            end
            A = A/(obj.detOutPlane);
            obj.psTensor.setValue(A);
        end
        
        function v = createVector(obj,d)
            v = zeros(d,1);
            if obj.isSymmbolic()
               v = sym(v); 
            end
        end
        
        function A = createSquareMatrix(obj,d)
            A = zeros(d,d);
            if obj.isSymmbolic()
                A = sym(A);
            end
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

