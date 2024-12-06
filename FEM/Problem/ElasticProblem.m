classdef ElasticProblem < handle
    
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

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFun();
            obj.createBCApplier();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacement();
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
            fun = {obj.uFun, obj.strainFun.project('P1',obj.mesh), ...
                obj.stressFun.project('P1',obj.mesh)};
            funNames = {'displacement', 'strain', 'stress'};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.scale       = cParams.scale;
            obj.mesh        = cParams.mesh;
            obj.solverType  = cParams.solverType;
            obj.solverMode  = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.uFun.ndimf;
            d.nnodes = size(obj.uFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end

        function createSolver(obj)
            sS.type      = obj.solverCase;
            solver       = Solver.create(sS);
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.bcApplier;
            obj.problemSolver    = ProblemSolver(s);
        end

        function computeStiffnessMatrix(obj)           
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.test     = obj.uFun;
            s.trial    = obj.uFun;
            s.material = obj.material;
            s.quadratureOrder = 2;
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
