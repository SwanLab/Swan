classdef StructuralTrussProblem < handle
    
    properties (Access = public)
        stress
        displacement
        adjoint
        stiffnessDerivative
    end

    properties (Access = private) % In
        coord
        connec
        material
        neumann, dirichlet
        barInterp, designVar
    end

    properties (Access = private) % Calc
        solver
        LHS, RHS
    end

    methods (Access = public)

        function obj = StructuralTrussProblem(cParams)
            obj.init(cParams);
            obj.createSolver();
        end

        function solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            obj.computeDisplacements();
        end

        function computeLHSDerivative(obj)
            s.coord    = obj.coord;
            s.connec   = obj.connec;
            s.material = obj.material;
%             s.
            lhs = LHSintegrator_StiffnessBeam(s);
            obj.stiffnessDerivative = lhs.computeDerivative();
        end

        function solveDisplacementAdjoint(obj,F)
            nNods = size(obj.coord, 1);
            nDofs = nNods * 3;
            fixedDofs = obj.convertDataToDof();
            fixedVals = fixedDofs(:,2);
            dofs = 1:nDofs;
            freeDofs = setdiff(dofs,fixedDofs(:,1));
            K   = obj.LHS(freeDofs, freeDofs);
            KRL = obj.LHS(freeDofs, fixedDofs(:,1));
            F = F - KRL * fixedVals;
            p = obj.solver.solve(K,F);
            obj.adjoint = p;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.coord     = cParams.coord;
            obj.connec    = cParams.connec;
            obj.neumann   = cParams.neumann;
            obj.dirichlet = cParams.dirichlet;
            obj.material  = cParams.material;
            obj.barInterp = cParams.interp;
            obj.designVar = cParams.designVar;
        end

        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeLHS(obj)
            % To optimize by calculating lengths just once
            % and just multiply the LHS by E*A/L
            s.coord    = obj.coord;
            s.connec   = obj.connec;
            s.material = obj.material;
            lhs = LHSintegrator_StiffnessBeam(s);
            obj.LHS = lhs.compute();
        end

        function computeRHS(obj)
            nNods = size(obj.coord, 1);
            nDofs = nNods * 3;
            s.neumann = obj.neumann;
            s.ndofs = nDofs;
            rhs = RHSintegrator_ISCSO(s);
            obj.RHS = rhs.compute();
        end

        function computeDisplacements(obj)
            nNods = size(obj.coord, 1);
            nDofs = nNods * 3;
            fixedDofs = obj.convertDataToDof();
            fixedVals = fixedDofs(:,2);
            dofs = 1:nDofs;
            freeDofs = setdiff(dofs,fixedDofs(:,1));
            K   = obj.LHS(freeDofs, freeDofs);
            KRL = obj.LHS(freeDofs, fixedDofs(:,1));
            F = obj.RHS(freeDofs) - KRL * fixedVals;
            u = obj.solver.solve(K,F);
            obj.displacement = u;
        end

        function Fdata = convertDataToDof(obj)
            Fdata = [];
            for i = 1:height(obj.dirichlet)
                data = obj.dirichlet(i,:);
                inod = data(1);
                iunk = data(2);
                idof = obj.nod2dof(inod,iunk);
                Fdata = [Fdata; idof, data(3)];
            end

        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = 3;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end
