classdef BoundaryConditionsPeriodic < BoundaryConditions

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
                    fullLHS  = obj.fullToReducedMatrix(K);
                    BCMatrix = fullLHS;
                    nConst   = 0;
            end
        end

        function fullRHS = computeBoundaryCondRHS(obj, cParams)
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
            perDOFslave = obj.periodic_constrained;
            obj.sizePer = size(perDOFslave, 1);
            if isprop(obj, 'vstrain')
                s.dirDOFs        = obj.dirichlet;
                s.ndofs          = obj.ndofs;
                DirComputer      = MacroDirichletComputer(s);   %a nivel de metodo igual micro que macro
                [CtDir, sizeDir] = DirComputer.computeDirCond();
                perDOFslave      = obj.periodic_constrained;
                perDOFmaster     = obj.periodic_free;
                obj.nConstraints = sizeDir + obj.sizePer; 
                Ct               = zeros(obj.nConstraints, obj.ndofs);
                for i = 1:obj.sizePer
                    masterNode                    = perDOFmaster(i);
                    slaveNode                     = perDOFslave(i);
                    Ct(i, [masterNode slaveNode]) = [1 -1];
                end
                Ct(obj.sizePer+1:end, :) = CtDir;
            else
                s.dirDOFs        = obj.dirichlet;
                s.ndofs          = obj.ndofs;
                DirComputer      = DirichletComputer(s);
                [CtDir, sizeDir] = DirComputer.computeDirCond();
                obj.nConstraints = sizeDir;
                Ct = CtDir;
            end
    
         end

         function fullRHS = createGeneralVector(obj, cParams) %Disp
                    
            if isprop(obj, 'vstrain')
                nPerDofs  = size(obj.periodic_constrained, 1);
                perVector = zeros(nPerDofs, 1);
                uD        = obj.dirichlet_values;
                fullRHS   = [cParams.RHS; perVector; uD]; 
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