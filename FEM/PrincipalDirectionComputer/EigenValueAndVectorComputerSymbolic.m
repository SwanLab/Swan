classdef EigenValueAndVectorComputerSymbolic < EigenValueAndVectorComputer
    
    properties (Access = private)
        eigenVectors
        eigenValues
    end
    
    methods (Access = public)
        
        function obj = EigenValueAndVectorComputerSymbolic(cParams)
            obj.init(cParams)
            obj.obtainEigenValuesAndVectors();
            obj.normalizeEigenVectors();
            obj.transformEigenVectorsInNormalBase();
            obj.transformSymEigenVectorsInFunctions();
            obj.transformSymEigenValuesInFunctions();
        end
        
    end
    
    methods (Access = private)
        
        function obtainEigenValuesAndVectors(obj)
             S = sym('S', [obj.ndim,obj.ndim],'real');
             S = triu(S,0) + triu(S,1).';
            [vS,dS] = eig(S);
            obj.eigenVectors = vS;
            obj.eigenValues = dS;
        end
        
        function normalizeEigenVectors(obj)
            for i = 1:obj.ndim
                e = obj.eigenVectors(:,i);
                obj.eigenVectors(:,i) = e/norm(e);
            end
        end
        
        function transformEigenVectorsInNormalBase(obj)
            %deter = det(obj.eigenVectors);
            %            obj.eigenVectors(:,2) = obj.eigenVectors(:,2)*deter;
            if size(obj.eigenVectors,1) == 2
                obj.eigenVectors(:,2) = [-obj.eigenVectors(2,1),obj.eigenVectors(1,1)];
            end
        end
        
        function transformSymEigenVectorsInFunctions(obj)
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    e = obj.eigenVectors(i,j);
                    obj.eigenVectorFunction{i,j} = matlabFunction(e);
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
    
end