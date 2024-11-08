classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        tau = 0.1
        tolerance = 1e-9

        type
        scale
        solverType
        solverMode
        solverCase
        quadOrder

        ElasticFun
        ExternalWorkFun

        u
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function results = compute(obj)
            bc = obj.boundaryConditions;
            uFun = LagrangianFunction.create(obj.mesh, 2, 'P1');
            uFun.fValues = obj.updateInitialDisplacement(bc,uFun);


            s.material = obj.material;
            s.u        = uFun;
            s.mesh     = obj.mesh;


            obj.ElasticFun = shFunc_Elastic(s);
            obj.ExternalWorkFun = shFunc_ExternalWork2(s);
            obj.u = uFun;
            
            errorU = 1;
            uOld = uFun.fValues;
            while (errorU >= obj.tolerance)

                LHS = obj.computeKInternal(obj.u);

                Fext = obj.computeFExternal(obj.u);
                Fint = obj.computeFInternal(obj.u);
                F = Fext+Fint;
                
                results = obj.computeU(F,obj.u,bc,LHS);
                 
                errorU = max(max(abs(uOld-results)));
                
                uOld = results;
                obj.u.fValues = results;
             end
    
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 2;
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material = cParams.material;
            obj.type      = cParams.type;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.solverCase = cParams.solverCase;
            obj.scale = cParams.scale;
        end

        function K = computeKInternal(obj,u)%,dispFun)
            K = obj.ElasticFun.computeHessian(obj.quadOrder,u);
        end

        function F = computeFInternal(obj,u)%,displacementFun,K)
            Ftry = obj.ElasticFun.computeJacobian(obj.quadOrder,u);
            F = Ftry;
        end


        function F = computeFExternal(obj,u)
            fExt = obj.boundaryConditions.pointloadFun;
            Ftry = obj.ExternalWorkFun.computeGradient(u,fExt,obj.quadOrder);
            F = Ftry;
        end

        function dim = getFunDims(obj,displacementFun)
            d.ndimf  = displacementFun.ndimf;
            d.nnodes = size(displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function uOut = computeU(obj,RHS,uIn,bc,LHS)


            % problemSolver = obj.createSolver();
            % 
            % t.stiffness = K;
            % t.forces = F;
            % u = problemSolver.solve(t);
            
            RHS = RHS(bc.free_dofs);
            LHS = LHS(bc.free_dofs,bc.free_dofs);

            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            uOutVec = uInVec;

            uInFree = uInVec(bc.free_dofs);
            uOutFree = obj.updateWithGradient(RHS,uInFree,LHS);
            uOutVec(bc.free_dofs) = uOutFree;
            uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';


        end
        % FOR INTERNAL MATERIA CREATION
        % function createMaterial (obj)
        %     s.type = 'ISOTROPIC';
        %     j = Material.create()
        %
        % end

        function problemSolver = createSolver(obj)
            sS.type      = obj.solverCase;
            solver       = Solver.create(sS);
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.createBCApplier();
            problemSolver    = ProblemSolver(s);
        end

        function BC = createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            BC = BCApplier(s);
        end


         function u = updateInitialDisplacement(obj,bc,uOld)
            restrictedDofs = bc.dirichlet_dofs;
            if isempty(restrictedDofs)
                u = uOld;
            else
                dirich = bc.dirichletFun;
                uVec = reshape(uOld.fValues',[uOld.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

                uVec(restrictedDofs) = dirichVec(restrictedDofs);
                u = reshape(uVec,[flip(size(uOld.fValues))])';
            end
         end

         function xNew = updateWithGradient(obj,RHS,x,LHS)
            % deltaX = -obj.tau.*RHS;
            deltaX = -LHS\RHS;
            xNew = x + deltaX; 
        end

    end
end