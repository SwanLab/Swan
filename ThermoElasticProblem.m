classdef ThermoElasticProblem < handle

    properties (Access = public)
        uFun
        strainFun
        stressFun
        forces

        tFun            
        tforces         
        T0
        alpha
    end

    properties (Access = private)
        quadrature
        boundaryConditionsElastic, bcApplierElastic

        stiffness
        solverType, solverMode, solverCase
        scale

        strain, stress

        problemSolverElastic

        kappa                             
        boundaryConditionsThermal, bcApplierThermal                               
        test                                        
        trial                                   
        source                                  
        problemSolverThermal                    
        Tstiffness                              
    end

    properties (Access = protected)
        mesh
        material
    end

    methods (Access = public)

        function obj = ThermoElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFun();
            obj.createBCApplierElastic();
            obj.createSolverElastic();
            
            obj.createTemperatureFun();           
            obj.createBCApplierThermal();         
            obj.createSolverThermal();            
        end

        function solve(obj)
            obj.computeThermalStiffnessMatrix();              
            obj.computeThermalForces();                      
            obj.computeTemperature();                         

            obj.computeStiffnessMatrix();                     
            obj.computeForces();                              
            obj.computeDisplacement();   

            obj.computeStrain();
            obj.computeStress();
        end

        function updateMaterial(obj, mat, kappa)
            obj.material = mat;
            obj.kappa = kappa;  
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

        function [tFun, uFun] = getTemperatureAndDisplacement(obj)
            tFun = obj.tFun;
            uFun = obj.uFun;
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.solverType  = cParams.solverType;
            obj.solverMode  = cParams.solverMode;
            obj.solverCase  = cParams.solverCase;

            % Elastic
            obj.scale       = cParams.scale;
            obj.material    = cParams.material;
            obj.boundaryConditionsElastic = cParams.boundaryConditionsElastic;

            % Thermal
            obj.alpha = cParams.alpha;
            obj.source       = cParams.source;                                    
            obj.T0           = cParams.T0;
            obj.boundaryConditionsThermal = cParams.boundaryConditionsThermal;    
            obj.test  = LagrangianFunction.create(obj.mesh,1,'P1');               
            obj.trial = LagrangianFunction.create(obj.mesh,1,'P1');               
        end

        
  %THERMAL PROBLEM (use thermalProblem problem)

        function createTemperatureFun(obj)
            obj.tFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function createBCApplierThermal(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditionsThermal;
            bcT = BCApplier(s);
            obj.bcApplierThermal = bcT;
        end

        function createSolverThermal(obj)
            s.solverType         = obj.solverType;
            s.solverMode         = obj.solverMode;
            s.solver             = obj.solverCase;
            s.boundaryConditions = obj.boundaryConditionsThermal;
            s.BCApplier          = obj.bcApplierThermal;
            obj.problemSolverThermal    = ProblemSolver(s);
        end

        function computeThermalStiffnessMatrix(obj)
            obj.Tstiffness=IntegrateLHS(@(u,v) obj.kappa.*DP(Grad(u),Grad(v)),obj.test,obj.trial,obj.mesh,'Domain',2);
        end

        function computeThermalForces(obj)
            obj.tforces = IntegrateRHS(@(v) DP(obj.source,v), obj.test, obj.mesh,'Domain',2);
        end

        function t = computeTemperature(obj)
            s.stiffness = obj.Tstiffness;
            s.forces    = obj.tforces;
            [t,~]       = obj.problemSolverThermal.solve(s);           
            obj.tFun.setFValues(t);
        end



%THERMO_ELASTIC PROBLEM

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createBCApplierElastic(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditionsElastic;
            bc = BCApplier(s);
            obj.bcApplierElastic = bc;
        end

        function createSolverElastic(obj)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solverCase;
            s.boundaryConditions = obj.boundaryConditionsElastic;
            s.BCApplier          = obj.bcApplierElastic;
            obj.problemSolverElastic    = ProblemSolver(s);
        end

        function computeStiffnessMatrix(obj)
            C     = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            obj.stiffness = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function computeForces(obj)
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
            obj.forces = rhs;

           %COUPLING TERM
            C     = obj.material;
            I = ConstantFunction.create(eye(2),obj.mesh);
            beta = obj.alpha.*DDP(C,I);       
            f = @(v) DDP(beta.*(obj.tFun - obj.T0),SymGrad(v));
            rhs_coupling = IntegrateRHS(f,obj.uFun,obj.mesh,'Domain',2);    
            rhs = rhs + rhs_coupling;
            obj.forces = rhs;
        end


        function u = computeDisplacement(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.problemSolverElastic.solve(s);
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
