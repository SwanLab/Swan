classdef PrincipalDirectionComputer < handle
   
   properties (Access = public)
       direction 
   end
   
   properties (Access = protected)
       directionFunction       
       ndim
       eigenVectors
   end    
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           f = PrincipalDirectionComputerFactory();
           obj = f.create(cParams);           
       end       
       
   end
   
   methods (Access = protected)
               
        function init(obj)
            obj.obtainEigenVectors();
            obj.normalizeEigenVectors();
            obj.transformSymEigenVectorsInFunctions();
        end       
               
        function normalizeEigenVectors(obj)
            for i = 1:obj.ndim
                e = obj.eigenVectors(:,i);
                obj.eigenVectors(:,i) = e/norm(e);
            end
        end
        
        function transformSymEigenVectorsInFunctions(obj)
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    e = obj.eigenVectors(i,j);
                    obj.directionFunction{i,j} = matlabFunction(e);
                end
            end
        end       
       
   end
   
   methods (Access = public, Abstract)
      compute(obj) 
   end   
    
   methods (Access = protected, Abstract)
      obtainEigenVectors(obj)
   end
end