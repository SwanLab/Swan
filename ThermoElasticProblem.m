classdef ThermoElasticProblem < handle

    properties (Access = public)
        uFun
        strainFun
        stressFun
        forces
    end

    properties (Access = private)
        quadrature
        boundaryConditions, bcApplier

        stiffness
        solverType, solverMode, solverCase
        scale

        strain, stress

        problemSolver
    end

    properties (Access = protected)
        mesh
        material
    end

    methods (Access = public)

        function obj = ThermoElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFun();
            obj.createBCApplier();
            obj.createSolver();
        end

        function solve(obj)
            % thermal problem
            % LHS thermal
            % RHS thermal
            % computeTemperature - obj.temperature 

            % for the thermo-elastic
            obj.computeStiffnessMatrix();   %LHS
            obj.computeForces();            %RHS - you need the temperature! 
            obj.computeDisplacement();      %Solve PDE
            obj.computeStrain();
            obj.computeStress();
        end

        function updateMaterial(obj, mat)
            obj.material = mat;
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
            fun = {obj.uFun, obj.strainFun.project('P1'), ...
                obj.stressFun.project('P1')};
            funNames = {'displacement', 'strain', 'stress'};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.scale       = cParams.scale;
            obj.solverType  = cParams.solverType;
            obj.solverMode  = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;

            % Temperature as a fixed function
%             T = LagrangianFunction.create(obj.mesh,1,'P1');
%             fValues = ones(Q.nDofs,1);
%             T.setFValues(fValues);
%             obj.temperature      = T;  
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end

        function createSolver(obj)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solverCase;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.bcApplier;
            obj.problemSolver    = ProblemSolver(s);
        end

        function computeStiffnessMatrix(obj)
            C     = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            obj.stiffness = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function computeForces(obj)  %this has to be changed
            bc  = obj.boundaryConditions;
            t   = bc.tractionFun;
            rhs = zeros(obj.uFun.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(obj.uFun);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditions;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -obj.stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(obj.uFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            %% - coupling term
            % f = @(v) beta*delta obj.temperature * div(v) 
            % rhs_coupling = IntegrateRHS(f,obj.testFun,obj.mesh,'Domain',2);    
            % rhs = rhs + rhs_coupling
            obj.forces = rhs;
        end


        function u = computeDisplacement(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.problemSolver.solve(s);
            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.uFun.setFValues(uSplit);
        end

        function computeStrain(obj)
            quad = Quadrature.create(obj.mesh, 2);
            xV = quad.posgp;
            obj.strainFun = SymGrad(obj.uFun);
            %             strFun = strFun.obtainVoigtFormat();
            obj.strain = obj.strainFun.evaluate(xV);
        end

        function computeStress(obj)
            quad = Quadrature.create(obj.mesh, 2);
            xV = quad.posgp;
            obj.stressFun = DDP(obj.material, obj.strainFun);
            obj.stress = obj.stressFun.evaluate(xV);
        end

    end

end
