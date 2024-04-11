classdef ElasticProblem < handle
    
    properties (Access = public)
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        quadrature
        boundaryConditions, BCApplier

        stiffness
        forces
        solver, solverType, solverMode
        scale
        
        strain, stress

        solverTol, solverParams, solverCase
        internalForcesComputer, testFunction
    end

    properties (Access = protected)
        mesh 
        material  
        displacementFun
    end

    methods (Access = public)

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDisplacementFun();
            obj.createBCApplier();
            obj.prepareInternalForcesFunction();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacement();
            % obj.computeStrain();
            % obj.computeStress();
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
            fun = {obj.uFun, obj.strainFun.project('P1',obj.mesh), ...
                obj.stressFun.project('P1',obj.mesh)};
            funNames = {'displacement', 'strain', 'stress'};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.material     = cParams.material;
            obj.scale        = cParams.scale;
            obj.mesh         = cParams.mesh;
            obj.solverType   = cParams.solverType;
            obj.solverMode   = cParams.solverMode;
            obj.solverTol    = cParams.solverTol;
            obj.solverParams = cParams.solverParams;
            obj.solverCase   = cParams.solverCase;
            obj.boundaryConditions = cParams.boundaryConditions;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createDisplacementFun(obj)
            obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.BCApplier = bc;
        end

        function createSolver(obj)
            s.tol          = obj.solverTol;
            s.solverParams = obj.solverParams;
            s.type         = obj.solverCase;
            obj.solver     = Solver.create(s);
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacementFun;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.stiffness = lhs.compute();
        end

        function computeForces(obj)
            s.type     = 'Elastic';
            s.scale    = 'MACRO';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(obj.solverType,'REDUCED')
                R = RHSint.computeReactions(obj.stiffness);
                obj.forces = rhs+R;
            else
                obj.forces = rhs;
            end
        end

        function u = computeDisplacement(obj)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solver;
            s.stiffness = obj.stiffness;
            s.forces = obj.forces;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier = obj.BCApplier;
            s.internalForces = @obj.computeInternalForces;
            pb = ProblemSolver(s);
            % - TO RECOVER
            [u,L] = pb.solve();
            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order   = 'P1';
            uFeFun = LagrangianFunction(z);
            obj.uFun = uFeFun;
            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function [u, L] = cleanupSolution(obj,sol)
            % bcapp = obj.BCApplier;
            bcs   = obj.boundaryConditions;
            % hasPeriodic = ~isequal(bcs.periodic_leader, []);
            switch true
                case strcmp(obj.solverType, 'MONOLITHIC')
                    nDisp = size(obj.stiffness,1);
                    u = sol(1:nDisp, :);
                    L = -sol( (nDisp+1):end, : );
                case strcmp(obj.solverType, 'REDUCED') && strcmp(obj.solverMode, 'DISP')
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    u = zeros(size(obj.stiffness,1), 1);
                    u(free_dofs) = sol;
                    u(bcs.dirichlet_dofs) = bcs.dirichlet_vals;
                    L = [];
                case strcmp(obj.solverType, 'REDUCED') && strcmp(obj.solverMode, 'FLUC')
                    lead = bcs.periodic_leader;
                    fllw = bcs.periodic_follower;
                    drch = bcs.dirichlet_dofs;
                    dofs = 1:size(obj.stiffness);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    u = zeros(length(dofs),1);
                    u(free) = sol(1:1:size(free,2));
                    u(lead) = sol(size(free,2)+1:1:size(sol,1));
                    u(fllw) = u(lead);
                    u(drch) = bcs.dirichlet_vals;
                    L = [];
                otherwise
                    u = [];
                    L = [];
            end
        end

        function prepareInternalForcesFunction(obj)
            s.type                     = 'ShapeSymmetricDerivative';
            s.scale                    = 'MACRO';
            s.mesh                     = obj.mesh;
            s.quadratureOrder          = obj.quadrature.order;
            s.globalConnec             = obj.mesh.connec;
            obj.internalForcesComputer = RHSintegrator.create(s);
            obj.testFunction           = LagrangianFunction.create(s.mesh,s.mesh.ndim,'P1');
        end

        function Fi = computeInternalForces(obj,x)
            [u,~]             = obj.cleanupSolution(x);
            uSplit            = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
            strainF           = SymGrad(obj.displacementFun);
            stressF           = DDP(obj.material, strainF);
            Fi                = obj.internalForcesComputer.compute(stressF,obj.testFunction);
            dofs              = 1:numel(Fi);
            free_dofs         = setdiff(dofs, obj.boundaryConditions.dirichlet_dofs);
            Fi                = Fi(free_dofs);
        end

        function computeStrain(obj)
            xV = obj.quadrature.posgp;
            obj.strainFun = SymGrad(obj.displacementFun);
%             strFun = strFun.obtainVoigtFormat();
            obj.strain = obj.strainFun.evaluate(xV);
        end

        function computeStress(obj)
            xV = obj.quadrature.posgp;
            obj.stressFun = DDP(obj.material, obj.strainFun);
            obj.stressFun.ndimf = obj.strainFun.ndimf;
            obj.stress = obj.stressFun.evaluate(xV);
        end

    end

end
