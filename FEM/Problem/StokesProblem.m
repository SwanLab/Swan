classdef StokesProblem < handle

    properties (Access = public)
        velocityFun
        pressureFun
    end
    
    properties (Access = private)
        mesh
        material
        solver

        state
        dtime
        finalTime
        inputBC
        boundaryConditions

        LHS, LHSintegrator
        RHS
    end

    methods (Access = public)

        function obj = StokesProblem(cParams)
            obj.init(cParams);
            obj.createVelocity();
            obj.createPressure();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.computeLHS();
            obj.computeRHS();
        end
        
        function computeVariables(obj)
            tol = 1e-6;
            bc  = obj.boundaryConditions;
            free_dof = [length(bc.freeFields{1}), length(bc.freeFields{2})];
            total_free_dof = sum(free_dof);
            LHSr = bc.fullToReducedMatrix(obj.LHS);
            RHSr = bc.fullToReducedVector(obj.RHS);

            switch obj.state
                case 'Steady'
                    x = obj.solver.solve(LHSr, RHSr);

                case 'Transient'
                    RHS0 = RHSr;
                    x_n(:,1) = zeros(total_free_dof,1);
                    x0 = zeros(total_free_dof,1);
                    
                    for istep = 2: obj.finalTime/obj.dtime
                        r = LHSr*x0 - RHSr;
                        while dot(r,r) > tol
                            inc_x = obj.solver.solve(LHSr, -r);
                            x = x0 + inc_x;
                            Fint = LHSr*x;
                            Fext = RHSr;
                            r = Fint - Fext;
                            x0 = x;
                        end
                        x_n(:,istep) = x;
                        RHSr = obj.updateRHS(RHS0, x0(1:free_dof(1)));
                    end
                    x = x_n;
            end
            fullx = obj.boundaryConditions.reducedToFullVector(x);
            vars = obj.separateVariables(fullx);
            obj.velocityFun.setFValues(obj.splitVelocity(vars.u));
            obj.pressureFun.setFValues(vars.p(:,end));
        end
       
        function print(obj, filename, software)
            if nargin == 2; software = 'Paraview'; end
            [fun, funNames] = obj.getFunsToPlot();
            a.mesh     = obj.mesh;
            a.filename = filename;
            a.fun      = fun;
            a.funNames = funNames;
            a.type     = software;
            pst = FunctionPrinter.create(a);
            pst.print();
        end

        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.velocityFun, obj.pressureFun};
            funNames = {'velocity', 'pressure'};
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.state      = cParams.state;
            obj.dtime      = cParams.dtime;
            obj.finalTime  = cParams.finalTime;
            obj.mesh       = cParams.mesh;
            obj.material   = cParams.material;
            inBC.pressure  = cParams.bc.pressure;
            inBC.velocity  = cParams.bc.velocity;
            inBC.pointload = [];
            inBC.velocityBC    = cParams.bc.velocityBC;
            inBC.forcesFormula = cParams.bc.forcesFormula;
            obj.inputBC    = inBC;
        end

        function createVelocity(obj)
            obj.velocityFun = LagrangianFunction.create(obj.mesh, 2, 'P2');
        end

        function createPressure(obj)
            obj.pressureFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function createBoundaryConditions(obj)
            vel = obj.velocityFun;
            prs = obj.pressureFun;

            borderCond = @(x) x(:,1) == 0 | x(:,1) == 1 | x(:,2) == 0 | x(:,2) == 1;
            dirDofs = obj.velocityFun.getDofsFromCondition(borderCond);
            dirich = obj.adaptForDirichletConditions(dirDofs);

            bcV.dirichlet = dirich;
            bcV.pointload = [];
            bcV.ndimf     = vel.ndimf;
            bcV.ndofs     = obj.velocityFun.nDofs;
            bcP.dirichlet = obj.inputBC.pressure;
            bcP.pointload = [];
            bcP.ndimf     = prs.ndimf;
            bcP.ndofs     = obj.pressureFun.nDofs;
            ndofs = obj.velocityFun.nDofs + obj.pressureFun.nDofs;
            s.dim   = [];
            s.scale = 'MACRO';
            s.bc    = {bcV, bcP};
            s.ndofs = ndofs; % Stokes
            bc = BoundaryConditionsStokes(s);
            obj.boundaryConditions = bc;
        end

        % Should be removed
        function dirich = adaptForDirichletConditions(obj,dirDofs)
            nodes = 1 + (dirDofs(2:2:end)-2)/obj.velocityFun.ndimf;
            nodes2 = repmat(nodes, [1 2]);
            iNod = sort(nodes2(:));
            mat12 = repmat([1;2], [length(iNod)/2 1]);
            valmat = zeros(length(iNod),1);
            dirich = [iNod mat12 valmat]; 
        end

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end
        
        function LHS = computeLHS(obj)
            s.type          = 'Stokes';
            s.dt            = obj.dtime;
            s.mesh          = obj.mesh;
            s.material      = obj.material;
            s.velocityFun = obj.velocityFun;
            s.pressureFun = obj.pressureFun;
            LHS_int = LHSintegrator.create(s);
            LHS = LHS_int.compute();
            obj.LHS = LHS;
            obj.LHSintegrator = LHS_int;
        end
        
        function RHS = computeRHS(obj)
            s.type          = 'Stokes';
            s.mesh          = obj.mesh;
            s.velocityFun   = obj.velocityFun;
            s.pressureFun   = obj.pressureFun;
            s.forcesFormula = obj.inputBC.forcesFormula;
            RHSint = RHSintegrator.create(s);
            F = RHSint.integrate();
            dirichlet = obj.boundaryConditions.dirichlet;
            uD = obj.boundaryConditions.dirichlet_values;
            R  = -obj.LHS(:,dirichlet)*uD;
            RHS = F + R;
            obj.RHS = RHS;
        end
        
        function variable = separateVariables(obj,x)
            ndofsV = obj.velocityFun.nDofs;
            variable.u = x(1:ndofsV,:);
            variable.p = x(ndofsV+1:end,:);
        end

        function RHS = updateRHS(obj, RHS0, x_n)
            RHS = zeros(size(RHS0));
            freeV = obj.boundaryConditions.freeFields{1};
            lenFreeV = length(freeV);
            M = obj.LHSintegrator.M;
            Mred = M(freeV,freeV);
            Mred_x_n = Mred*x_n;
            RHS(1:lenFreeV,1) = RHS0(1:lenFreeV,1) + Mred_x_n;
        end

        function uM = splitVelocity(obj, vel)
            u = vel;
            nu = obj.velocityFun.ndimf;
            nnode = round(length(u)/nu);
            nodes = 1:nnode;
            uM = zeros(nnode,nu);
            for idim = 1:nu
                dofs = nu*(nodes-1)+idim;
                uM(:,idim) = u(dofs, end);
            end
        end

    end

end
