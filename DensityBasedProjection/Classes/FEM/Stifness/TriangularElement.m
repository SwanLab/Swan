classdef TriangularElement < Element
    methods (Access = public)
        function computeStifnessMatrix(obj)
            A11     = [2  1 -1;  1  2  1; -1  1  2];
            A12     = [-1 -1  1;  1  0 -1; -1  1  0];
            B11     = [1 -1  0; -1  1  0;  0  0  0];
            B12     = [0  0  0;  0 -1  1;  0  1 -1];

            obj.stifnessMatrix = obj.t/(2*(1+obj.poissonCoefficient))*[A11 A12; A12' B11] + obj.t*obj.poissonCoefficient/(2*(1+obj.poissonCoefficient))*[B12 B12'; A12' A11];
        end
    end
end