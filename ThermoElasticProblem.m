classdef ThermoElasticProblem < handle

    properties (Access = public)
        tFun  %

        uFun
        strainFun
        stressFun
        forces
    end

    properties (Access = private)
        quadrature
        boundaryConditionsThermal
        boundaryConditionsElastic

        stiffness
        solverType, solverMode, solverCase

        strain, stress

        thermalProblemSolver
        thermoElasticProblemSolver

        conductivity  %
        test          %
        trial         %
        source        %
        temperature   %
    end

    properties (Access = protected)
        mesh
        materialElastic
    end

    methods (Access = public)

        function obj = ThermoElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFun();
            obj.createTemperatureFun();    %
            obj.createThermalSolver();
            obj.createThermoElasticSolver();
        end

        function solve(obj)
            % thermal problem
            % LHS thermal
            % RHS thermal
            % computeTemperature - obj.temperature 
            obj.computeThermalStiffnessMatrix(kappa); % LHS
            obj.computeThermalForces();          % RHS
            obj.computeTemperature();     % Solve PDE 

            % for the thermo-elastic
            obj.computeStiffnessMatrix();   %LHS
            obj.computeForces();            %RHS - you need the temperature! 
            obj.computeDisplacement();      %Solve PDE
            obj.computeStrain();
            obj.computeStress();
        end

        function updateMaterial(obj, mat)
            obj.materialElastic = mat;
        end

        % updateExternalForces... Ask Giovanna

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
            obj.mesh = cParams.mesh;
            obj.materialElastic = cParams.material;
            % Interpolator thermal problem
            obj.solverType  = cParams.solverType;
            obj.solverMode  = cParams.solverMode;
            obj.solverCase  = cParams.solverCase;
            obj.boundaryConditionsThermal = cParams.boundaryConditionsThermal;
            obj.boundaryConditionsElastic = cParams.boundaryConditionsElastic;

            obj.test  = LagrangianFunction.create(obj.mesh,1,'P1');   %
            obj.trial = LagrangianFunction.create(obj.mesh,1,'P1');  %

            % Temperature as a fixed function
            T = LagrangianFunction.create(obj.mesh,1,'P1');
            fValues = ones(T.nDofs,1);
            T.setFValues(fValues);
            obj.temperature      = T;  
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function bcAp = createBCApplier(obj,bc)
            s.mesh = obj.mesh;
            s.boundaryConditions = bc;
            bcAp = BCApplier(s);
        end

        function createThermalSolver(obj)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solverCase;
            s.boundaryConditions = obj.boundaryConditionsThermal;
            s.BCApplier          = obj.createBCApplier(obj.boundaryConditionsThermal);
            obj.thermalProblemSolver    = ProblemSolver(s);
        end

        function createThermoElasticSolver(obj)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solverCase;
            s.boundaryConditions = obj.boundaryConditionsElastic;
            s.BCApplier          = obj.createBCApplier(obj.boundaryConditionsElastic);
            obj.thermoElasticProblemSolver    = ProblemSolver(s);
        end

        function computeStiffnessMatrix(obj)
            C     = obj.materialElastic;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            obj.stiffness = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function computeForces(obj)  %this has to be changed
            bc  = obj.boundaryConditionsElastic;
            t   = bc.tractionFun;
            rhs = zeros(obj.uFun.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(obj.uFun);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditionsElastic;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -obj.stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(obj.uFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            % - coupling term
            f = @(v) beta*obj.temperature * div(v); 
            rhs_coupling = IntegrateRHS(f,obj.test,obj.mesh,'Domain',2);    
            rhs = rhs + rhs_coupling;
            obj.forces = rhs;
        end

        % THERMAL PROBLEM;

        function createTemperatureFun(obj)
            obj.tFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

         function computeThermalStiffnessMatrix(obj, kappa)
            obj.stiffness=IntegrateLHS(@(u,v) kappa.*DP(Grad(u),Grad(v)),obj.test,obj.trial,obj.mesh,'Domain',2);
        end

        function computeThermalForces(obj)
            obj.forces = IntegrateRHS(@(v) DP(obj.source,v), obj.test, obj.mesh,'Domain',2);
        end

        function u = computeTemperature(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.thermalProblemSolver.solve(s);           
            obj.tFun.setFValues(u);
        end



        % ELASTIC COMPUTATION
        function u = computeDisplacement(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.thermoElasticProblemSolver.solve(s);
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
            obj.stressFun = DDP(obj.materialElastic, obj.strainFun);
            obj.stress = obj.stressFun.evaluate(xV);
        end

    end

end
