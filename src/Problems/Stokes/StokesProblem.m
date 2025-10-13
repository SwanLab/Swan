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

            vF  = obj.velocityFun;
            v0F = LagrangianFunction.create(obj.mesh,vF.ndimf,vF.order);


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

                        freeV = obj.boundaryConditions.freeFields{1};

                        x0T = obj.boundaryConditions.reducedToFullVector(x0);
                        x0Tv = x0T(1:obj.velocityFun.nDofs);
                        v0F.setFValues(obj.splitVelocity(x0Tv));

                        Fext = IntegrateRHS(@(v) DP(v,v0F)./obj.dtime,vF,obj.mesh,3);            
                        RHSt = Fext(freeV);                        

                        RHSr = zeros(size(RHS0));
                        freeV = obj.boundaryConditions.freeFields{1};
                        lenFreeV = length(freeV);

                        RHSr(1:lenFreeV,1) = RHS0(1:lenFreeV,1) + RHSt;

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
            vF = obj.velocityFun;
            pF = obj.pressureFun;

            borderCond = @(x) x(:,1) == 0 | x(:,1) == 1 | x(:,2) == 0 | x(:,2) == 1;
            dirDofs = obj.velocityFun.getDofsFromCondition(borderCond);
            dirich = obj.adaptForDirichletConditions(dirDofs);

            bcV.dirichlet = dirich;
            bcV.pointload = [];
            bcV.ndimf     = vF.ndimf;
            bcV.ndofs     = obj.velocityFun.nDofs;
            bcP.dirichlet = obj.inputBC.pressure;
            bcP.pointload = [];
            bcP.ndimf     = pF.ndimf;
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
        
        function computeLHS(obj)
            vF = obj.velocityFun;
            uF = vF;
            pF = obj.pressureFun;
            M = IntegrateLHS(@(u,v) DP(v,u)./obj.dtime,vF,uF,obj.mesh,3);
            K = IntegrateLHS(@(u,v) DDP(Grad(v),Grad(u)),vF,uF,obj.mesh);            
            D = IntegrateLHS(@(p,v) DP(v,Grad(p)),vF,pF,obj.mesh);
            Z = sparse(size(D, 2),size(D, 2));
            obj.LHS = [K+M, D; D',Z];
        end
        
        function computeRHS(obj)
            f = obj.inputBC.forcesFormula;
            vF = obj.velocityFun;
            pF = obj.pressureFun;
            Fext = IntegrateRHS(@(v) DP(v,f),vF,obj.mesh,3);            
            F = [Fext; zeros(pF.nDofs,1)];



            dirichlet = obj.boundaryConditions.dirichlet;
            uD = obj.boundaryConditions.dirichlet_values;
            R  = -obj.LHS(:,dirichlet)*uD;



            obj.RHS = F + R;
        end
        
        function variable = separateVariables(obj,x)
            ndofsV = obj.velocityFun.nDofs;
            variable.u = x(1:ndofsV,:);
            variable.p = x(ndofsV+1:end,:);
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
