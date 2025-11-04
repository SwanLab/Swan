classdef Newton
    methods
        function [uOut, uOutVec] = computeDisplacement(obj, LHSfull, RHSfull, uIn, bc)
            [LHS, RHS] = obj.fullToReduced(LHSfull, RHSfull, bc);

            if ~isempty(LHS)
                uInVec = reshape(uIn.fValues', [uIn.nDofs 1]);
                uOutVec = uInVec;

                uInFree = uInVec(bc.free_dofs);
                uOutFree = obj.updateWithNewton(LHS, RHS, uInFree);

                uOutVec(bc.free_dofs) = uOutFree;
                uOut = reshape(uOutVec, [flip(size(uIn.fValues))])';
            else
                uOut = uIn.fValues;
                uOutVec = reshape(uIn.fValues', [uIn.nDofs 1]);
            end
        end
    end

    methods (Access = private)
        function [LHS, RHS] = fullToReduced(~, LHS, RHS, bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function xNew = updateWithNewton(~, LHS, RHS, x)
            deltaX = -LHS \ RHS;
            xNew = x + deltaX;
        end
    end
end
