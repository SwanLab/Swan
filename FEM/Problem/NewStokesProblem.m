classdef NewStokesProblem < handle

    properties (Access = public)
        variables
    end
    
    properties (Access = private)
        interp
        problemData
        geometry
        mesh
        dof
        element
        material
        solver
        fileName
    end

    methods (Access = public)

        function obj = NewStokesProblem(fileName)
            obj.fileName = fileName;
            obj.readProblemData(fileName);
            obj.createGeometry();
            obj.createInterpolation();
            obj.createDOF();
            obj.createBoundaryConditions();
            obj.createMaterial();
            obj.createElement();
            obj.createSolver();
        end
        
        function computeVariables(obj)
            nFields = numel(obj.interp);
            for ifield = 1:nFields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            transient = false;  % !! This should not be defined in here !!
            tol = 1e-6;         % !! This should not be defined in here !!
            if transient
                dt = 0.01;      % !! This should not be defined in here !!
                final_time = 1; % !! This should not be defined in here !!
                x = obj.solveTransientNonlinearProblem(free_dof,tol,dt,final_time);
            else
                x = obj.solveSteadyNonlinearProblem(free_dof,tol);
            end
            obj.variables = obj.element.computeVars(x);
        end

    end
    
    methods (Access = private)
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry    = Geometry.create(s);
            obj.geometry(2) = Geometry.create(s);
        end
        
        function createInterpolation(obj)
            interpU = 'QUADRATIC';
            interpP = 'LINEAR';
            obj.interp{1}=Interpolation.create(obj.mesh,interpU);
            obj.interp{2}=Interpolation.create(obj.mesh,interpP);
        end

        function createDOF(obj)
            obj.dof = DOF_Stokes(obj.fileName,obj.mesh,obj.geometry,obj.interp);
        end

        function createBoundaryConditions(obj)
        end

        function createMaterial(obj)
            cParams.nelem = obj.mesh.nelem;
            obj.material = Material_Stokes(cParams);
        end

        function createElement(obj)
            obj.element  = Element_Stokes(obj.geometry,obj.mesh,obj.material,obj.dof,obj.problemData,obj.interp);
        end

        function createSolver(obj)
            obj.solver   = Solver.create();
        end

        function sol = solveSteadyNonlinearProblem(obj,free_dof,tol)
            total_free_dof = sum(free_dof);
            dr = obj.element.computedr;
            x0 = zeros(total_free_dof,1);
            
            r = obj.element.computeResidual(x0,dr);
            x = x0;
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r
                r = obj.element.computeResidual(x,dr);
                x0 = x;
            end
            sol = x;
        end
        
        function sol = solveTransientNonlinearProblem(obj,free_dof,tol,dt,final_time)
            total_free_dof = sum(free_dof);
            x_n(:,1) = zeros(total_free_dof,1);
            x0 = zeros(total_free_dof,1);
            
            dr = obj.element.computedr(dt);
            
            for istep = 2: final_time/dt
                u_previous_step = x_n(1:free_dof(1),istep-1);
                
                r = obj.element.computeResidual(x0,dr,u_previous_step);
                while dot(r,r) > tol
                    inc_x = obj.solver.solve(dr,-r);
                    x = x0 + inc_x;
                    % Compute r
                    r = obj.element.computeResidual(x,dr,u_previous_step);
                    x0 = x;
                end
                x_n(:,istep) = x;
            end
            sol = x_n;
        end
        
        function readProblemData(obj,fileName)
            femReader = FemInputReader_GiD();
            s = femReader.read(fileName);
            
            obj.problemData.fileName = fileName;
            obj.problemData.scale = s.scale;
            obj.problemData.pdim  = s.pdim;
            obj.problemData.ptype = s.ptype;
            obj.problemData.nelem = s.mesh.nelem;
            obj.problemData.bc.dirichlet = s.dirichlet;
            obj.problemData.bc.pointload = s.pointload;
            obj.mesh = s.mesh;
        end

    end

end
