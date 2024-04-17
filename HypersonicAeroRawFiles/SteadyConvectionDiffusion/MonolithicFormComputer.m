classdef MonolithicFormComputer < handle

    properties (Access = private)
        stiffnessMatrix
        extFluxesVector
        Cmatrix
        bc
    end

    methods (Access = public)
        function obj = MonolithicFormComputer(cParams)
            obj.init(cParams);
        end

        function [Ktot,ftot] = compute(obj)
            b    = obj.bc.dirichlet_vals;
            K    = obj.stiffnessMatrix;
            f    = obj.extFluxesVector;
            A    = obj.Cmatrix;
            Ktot = [K A';A zeros(2,2)];
            ftot = [f;b]; % There is an issue in P2 case !!
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.stiffnessMatrix = cParams.K;
            obj.extFluxesVector = cParams.f;
            obj.bc              = cParams.bc;
            obj.Cmatrix = obj.computeMultipliersCoefficientsMatrix(cParams);
        end

        function A = computeMultipliersCoefficientsMatrix(obj,cParams)
            dirDofs = obj.bc.dirichlet_dofs;
            switch cParams.p
                case 1
                    numnp = cParams.nnodes;
                case 2
                    numnp = 2*cParams.nnodes - 1;
            end
            A            = zeros(2,numnp);
            A(1,dirDofs(1)) = 1;
            A(2,dirDofs(2)) = 1;
        end
    end
end