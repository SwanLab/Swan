classdef BoundaryConditionsLineal < BoundaryConditions

    properties (Access = private)
        nConstraints
        sizePer
    end

    methods (Access = public)
        function [BCMatrix, nConst] = computeBoundaryCondLHS(obj, K)
            switch obj.solType 
                case 'MONOLITIC'
                    Ct       = obj.createConstraintMatrix();
                    BCMatrix = Ct;
                    nConst   = obj.nConstraints;
                case 'REDUCED'
                    fullLHS  = obj.fullToReducedMatrix(K); % nope should be in solver
                    BCMatrix = fullLHS;
                    nConst   = 0;
            end
        end

        function fullRHS  = computeBoundaryCondRHS(obj, cParams)
            switch obj.solType 
                case 'MONOLITIC'
                    fullRHS     = obj.createGeneralVector(cParams);
                case 'REDUCED'
                    R           = obj.computeReactions(cParams.LHS);
                    CompleteRHS = cParams.RHS + R;
                    fullRHS     = obj.fullToReducedVector(CompleteRHS);
            end
        end
    end

    methods (Access = private)
        function Ct = createConstraintMatrix(obj)
            if ~isempty(obj.vstrain)
                perDOFslave = obj.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
                [CtDir, CtPerDir] = obj.conditionSetter.setDirichletLhs();
                perDOFmaster      = obj.periodic_free;
                perDOFslave       = obj.periodic_constrained;
                nEqperType        = obj.sizePer/4;
                obj.nConstraints  = size(CtDir, 1)+ size(CtPerDir, 1) ...
                                    + nEqperType*3; 
                Ct                = zeros(nEqperType*3+size(CtDir, 1)+ ...
                                    size(CtPerDir, 1), obj.ndofs);
                for i = 1:nEqperType
                    masterDOF                   = perDOFmaster(i);
                    slaveDOF                    = perDOFslave(i);
                    Ct(i, [masterDOF slaveDOF]) = [1 -1];
                end
                for i = nEqperType+1:2*nEqperType
                    masterDOF                   = perDOFmaster(i);
                    slaveDOF                    = perDOFslave(i);
                    Ct(i, [masterDOF slaveDOF]) = [1 -1];
                end
                for i = 2*nEqperType+1:3*nEqperType
                    masterDOF                            = perDOFmaster(i);
                    slaveDOF                             = perDOFslave(i);
                    positionCt                           = i-nEqperType;
                    Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
                end
                for i = 3*nEqperType+1:4*nEqperType
                    masterDOF                            = perDOFmaster(i);
                    slaveDOF                             = perDOFslave(i);
                    positionCt                           = i-nEqperType;
                    Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
                end
                Ct(nEqperType*3+1:end, :) = [CtPerDir; CtDir];
            else
                s.dirDOFs        = obj.dirichlet;
                s.ndofs          = obj.ndofs;
                DirComputer      = MacroDirichletComputer(s);
                [CtDir, sizeDir] = DirComputer.computeDirCond;
                obj.nConstraints = sizeDir;
                Ct = CtDir;
            end

        end

        function fullRHS = createGeneralVector(obj, cParams) 
            if ~isempty(obj.vstrain)
                nU           = size(cParams.LHS, 1);
                der1Rhs      = zeros(nU, 1);
                nEqperType   = obj.sizePer/4;
                perVector    = zeros(nEqperType*3, 1);
                for i = 1:nEqperType
                    perVector(i, 1) = -obj.vstrain(1);
                end
                for i = nEqperType+1:2*nEqperType
                    perVector(i, 1) = -obj.vstrain(3);
                end
                for i = 2*nEqperType+1:3*nEqperType
                    perVector(i, 1) = -obj.vstrain(2);
                end
                [RHSDir, RHSDirPer] = obj.conditionSetter.setDirichletRhs();
                fullRHS             = [der1Rhs; perVector; RHSDirPer; RHSDir];
            else
                uD        = obj.dirichlet_values;
                fullRHS   = [cParams.RHS; uD];
            end

        end

        function R = computeReactions(obj, K)
            dirich        = obj.dirichlet;
            dirichV       = obj.dirichlet_values;
            if ~isempty(dirich)
                R = -K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndofs(:)),1);
            end
        end

    end
end