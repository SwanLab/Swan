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
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacement();
            obj.computeStrain();
            obj.computeStress();
        end

        function plot(obj)
            s.dim          = obj.getFunDims();
            s.mesh         = obj.mesh;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.getFunDims();
        end

        function updateMaterial(obj, mat)
            obj.material = mat;
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

        function quad = getQuadrature(obj)
            quad  = obj.quadrature;
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
            fun = {obj.uFun, obj.strainFun, obj.stressFun};
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
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
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
            s.stiffness = obj.stiffness;
            s.forces = obj.forces;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier = obj.BCApplier;
            pb = ProblemSolver(s);
            [u,L] = pb.solve();
            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order   = 'P1';
            uFeFun = LagrangianFunction(z);
            obj.uFun = uFeFun;
            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrain(obj)
            xV = obj.quadrature.posgp;

            strFun = SymGrad(obj.displacementFun).evaluate(xV);
%             strFun = strFun.obtainVoigtFormat();
            obj.strainFun = strFun;
            obj.strain = strFun;
        end

        function computeStress(obj)
            strn(:,1,:,:) = obj.strain;
            Cv            = obj.material.evaluate(obj.quadrature.posgp);

            strs = pagemtimes(Cv,strn);
            strs = permute(strs, [1 3 4 2]);

            z.mesh       = obj.mesh;
            z.fValues    = strs;
            z.quadrature = obj.quadrature;
%             strFun       = FGaussDiscontinuousFunction(z);

%             obj.stress    = strFun;
%             obj.stressFun = strFun;
        end

    end

end
