classdef NewStokesProblem < handle

    properties (Access = public)
        variables
    end
    
    properties (Access = private)
        mesh
        material
        solver

        dim
        state
        dtime
        finalTime
        inputBC
        velocityField
        pressureField
        boundaryConditions

        LHS, LHSintegrator
        RHS
    end

    methods (Access = public)

        function obj = NewStokesProblem(cParams)
            obj.init(cParams);
            obj.createVelocityField();
            obj.createPressureField();
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
                    obj.variables = obj.separateVariables(x);

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
                    obj.variables = obj.separateVariables(x);
            end
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

        function createVelocityField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 2;
            s.interpolationOrder = 'QUADRATIC';
            s.scale              = 'MACRO';
            obj.velocityField = Field(s);
        end

        function createPressureField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            s.scale              = 'MACRO';
            obj.pressureField = Field(s);
        end

        function createBoundaryConditions(obj)
            vel = obj.velocityField;
            prs = obj.pressureField;
            newBC = vel.adaptBoundaryConditions(obj.inputBC);
            bcV.dirichlet = newBC.dirichlet;
            bcV.pointload = [];
            bcV.ndimf     = vel.dim.ndimf;
            bcV.ndofs     = vel.dim.ndofs;
            bcP.dirichlet = obj.inputBC.pressure;
            bcP.pointload = [];
            bcP.ndimf     = prs.dim.ndimf;
            bcP.ndofs     = prs.dim.ndofs;
            ndofs = vel.dim.ndofs + prs.dim.ndofs;
            s.dim   = [];
            s.scale = 'MACRO';
            s.bc    = {bcV, bcP};
            s.ndofs = ndofs; % Stokes
            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
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
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            LHS_int = LHSintegrator.create(s);
            LHS = LHS_int.compute();
            obj.LHS = LHS;
            obj.LHSintegrator = LHS_int;
        end
        
        function RHS = computeRHS(obj)
            s.type          = 'Stokes';
            s.mesh          = obj.mesh;
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            s.forcesFormula = obj.inputBC.forcesFormula;
            RHSint = RHSintegrator.create(s);
            F = RHSint.integrate();
            dirichlet = obj.boundaryConditions.dirichlet;
            uD = obj.boundaryConditions.dirichlet_values;
            R  = -obj.LHS(:,dirichlet)*uD;
            RHS = F + R;
            obj.RHS = RHS;
        end
        
        function variable = separateVariables(obj,x_free)
            x = obj.boundaryConditions.reducedToFullVector(x_free);
            ndofsV = obj.velocityField.dim.ndofs;
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

    end

end
