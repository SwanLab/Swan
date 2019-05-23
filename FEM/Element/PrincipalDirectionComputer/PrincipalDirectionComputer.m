classdef PrincipalDirectionComputer < handle
   
   properties (Access = public)
       direction 
       principalStress
   end
   
   properties (Access = protected)
       directionFunction       
       eigenValueFunction
       ndim
       eigenVectors
       eigenValues
       avarageTensor       
   end    
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           f = PrincipalDirectionComputerFactory();
           obj = f.create(cParams);           
       end       
       
   end
   
   methods (Access = protected)
               
        function init(obj)
            obj.obtainEigenValuesAndVectors();
            obj.normalizeEigenVectors();
            obj.transformSymEigenVectorsInFunctions();
            obj.transformSymEigenValuesInFunctions();            
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
        
        function transformSymEigenValuesInFunctions(obj)
            for i = 1:obj.ndim
                    e = obj.eigenValues(i,i);
                    obj.eigenValueFunction{i} = matlabFunction(e);
            end            
        end
       
   end
   
   methods (Access = public, Abstract)
      compute(obj) 
   end   
    
   methods (Access = protected, Abstract)
      obtainEigenValuesAndVectors(obj)
   end
end