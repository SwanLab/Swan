classdef TriangularElement < Element 
    methods (Access = public)
        function computeStifnessMatrix(obj)           
            A = [1/2, 1/2, 0; 0, 1/2, 1/2];
            B = [-1, 1, 0; -1, 0, 1];
            D = obj.t/(1-obj.poissonCoefficient^2)*[1, obj.poissonCoefficient, 0; obj.poissonCoefficient, 1, 0; 0, 0, (1-obj.poissonCoefficient)/2];
            T = abs(det([1, obj.nodes(1, :); 1, obj.nodes(2, :); 1, obj.nodes(3, :)]))/2;
            obj.stiffnessMatrix = T*D*(B'/(A*B))*B'/(A*B)*A;
        end
    end
end