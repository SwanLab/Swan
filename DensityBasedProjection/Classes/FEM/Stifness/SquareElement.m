classdef SquareElement < handle 
    properties (Access = public)
        stifnessMatrix
    end
    properties(Access = private)
        poissonCoefficient
        t
    end
    

    methods (Access = public)
        function obj = SquareElement(cParams)
            obj.inputData(cParams);
        end

        function computeStifnessMatrix(obj)
            A11     = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
            A12     = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
            B11     = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
            B12     = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
            obj.stifnessMatrix      = obj.t/(1-obj.poissonCoefficient^2)/24*([A11 A12;A12' A11]+obj.poissonCoefficient*[B11 B12;B12' B11]);
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.poissonCoefficient = cParams.poissonCoefficient; 
            obj.t = cParams.t; 
        end
    end
end