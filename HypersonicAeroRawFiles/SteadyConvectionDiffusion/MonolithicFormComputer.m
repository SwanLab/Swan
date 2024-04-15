classdef MonolithicFormComputer < handle

    properties (Access = private)
        stiffnessMatrix
        extFluxesVector
        Cmatrix
    end

    methods (Access = public)
        function obj = MonolithicFormComputer(cParams)
            obj.init(cParams);
        end

        function [Ktot,ftot] = compute(obj,problem)
            switch problem
                case 2
                    b = [0;0];
                otherwise
                    b = [0;1];
            end
            K    = obj.stiffnessMatrix;
            f    = obj.extFluxesVector;
            A    = obj.Cmatrix;
            Ktot = [K A';A zeros(2,2)];
            ftot = [f;b];
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.Cmatrix = obj.computeMultipliersCoefficientsMatrix(cParams);
            obj.stiffnessMatrix = cParams.K;
            obj.extFluxesVector = cParams.f;
        end
    end

    methods (Static,Access = private)
        function A = computeMultipliersCoefficientsMatrix(cParams)
            nnodes      = cParams.nnodes;
            A           = zeros(2,nnodes);
            A(1,1)      = 1;
            A(2,nnodes) = 1;
        end
    end
end