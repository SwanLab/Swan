classdef ElasticProblem < handle
    
    properties (Access = public)
        variables
        boundaryConditions
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        stiffness
        forces
        solver
        integratorBuilder
        geometry
        scale
        inputBC
        strain
        stress

        newBC
    end

    properties (Access = protected)
        solType
        solMode
        quadrature
        material
        btype
        vstrain
        mesh % For Homogenization
        interpolationType
        displacementFun
    end

    methods (Access = public)

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFun();
            obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.compDisp();
            obj.computeDisplacements();
            obj.computeStrain();
            obj.computeStress();
            obj.computePrincipalDirection();
        end

        function plot(obj)
            s.dim          = obj.getFunDims();
            s.mesh         = obj.mesh;
%             s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.getFunDims();
        end

        function setC(obj, C)
            obj.material.C = C;
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

        function quad = getQuadrature(obj)
            quad  = obj.quadrature;
        end
       
        function print(obj, filename, software)
            if nargin == 2; software = 'GiD'; end
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
            obj.inputBC     = cParams.bc;
            obj.newBC = cParams.newBC;

%             obj.btype       = cParams.builderType;
            obj.solMode = cParams.solMode;
            obj.solType = cParams.solType;                        
%             if isprop(cParams, 'interpolationType')

            obj.mesh        = cParams.mesh;
            if isfield(cParams, 'interpolationType')
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end
            obj.createQuadrature();
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createDisplacementFun(obj)
            obj.displacementFun = P1Function.create(obj.mesh, obj.mesh.ndim);
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBoundaryConditions(obj)
            dim = obj.getFunDims();
            bc = obj.inputBC;
            bc.ndimf = dim.ndimf;
            bc.ndofs = dim.ndofs;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = dim.ndofs;
            s.btype = obj.btype; % builder type
            s.solMode = obj.solMode;
            s.solType = obj.solType;
            if ~isempty(obj.vstrain)
                s.vstrain = obj.vstrain;
            end
            f = BoundaryConditionsFactory();
            bc = f.create(s);
%             s.ndofs = dim.ndofs;
%             bc = BoundaryConditions(s);
            bc.compute();
            %MOVER a computeDisplacements para que vstrain esté definido.
            obj.boundaryConditions = bc;
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
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;

            s.globalConnec = obj.mesh.connec;
            s.btype = obj.btype;

%             s.globalConnec = obj.displacementField.connec;
            s.globalConnec = obj.mesh.connec;


            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();

            obj.variables.fext = rhs;

            R = RHSint.computeReactions(obj.stiffness);
%             obj.variables.fext = rhs + R;

            obj.forces = rhs + R;
        end
        
        function u = compDisp(obj)
            s.type = 'MONOLITHIC';
            s.stiffness = obj.stiffness;
            s.forces = obj.forces;
            s.boundaryConditions = obj.newBC;
            pb = ProblemSolver(s);
            u = pb.solve();
            % u = 1;
            % u = ProblemSolver.solve(LHS,RHS, 'MONOLITHIC');
        end

        function u = computeDisplacements(obj)

%             bc = obj.boundaryConditions;
%             Kred = bc.fullToReducedMatrix(obj.LHS);
%             Fred = bc.fullToReducedVector(obj.RHS);
%             u = obj.solver.solve(Kred,Fred);
%             u = bc.reducedToFullVector(u);
%             obj.variables.d_u = u;


            bc            = obj.boundaryConditions;
%             builder       = bc.integratorBuilder;
%             s.LHS         = obj.stiffnessMatrix;
            s.LHS         = obj.stiffness;
            s.RHS         = obj.forces;
            s.bc          = bc;
            s.builderType = obj.btype;
            s.solver      = obj.solver;
            s.scale       = obj.scale;
            s.mesh        = obj.mesh;
            dim = obj.getFunDims();
            s.dim         = dim;
            
            % ConstraintSolver
            s.RHS         = obj.forces;
            s.bc          = obj.boundaryConditions;
            s.K           = obj.stiffness;
            s.solver      = obj.solver;
%             if ~isempty(obj.vstrain)
%                 s.vstrain = obj.vstrain;
%             end

            s.solType = obj.solType;
            
            %REFORMULATED ARCHITECTURE
            if strcmp(obj.solType, 'MONOLITIC') && strcmp(obj.solMode, 'DISP')  
                    bc.computeMonoliticMicroConditionDisp();
            end

%             f = ConstraintSolverFactory();
            bcSolver = ConstraintSolverFactory.create(s);
            [LHSMatrix, nConstraints] = bc.computeBoundaryCondLHS(s.LHS);
            RHSMatrix = bc.computeBoundaryCondRHS(s);
            [u, L] = bcSolver.solveSystem(LHSMatrix, RHSMatrix, nConstraints);
            obj.variables.d_u = u;
            obj.variables.LangMult = L;
            
            z.mesh   = obj.mesh;


            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFeFun = P1Function(z);
            obj.uFun = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrain(obj)
            strFun = obj.displacementFun.computeSymmetricGradient(obj.quadrature);
            strFun.applyVoigtNotation();
            perm = permute(strFun.fValues, [2 1 3]);
            obj.variables.strain = perm;
            obj.strainFun = strFun;
            obj.strain = strFun;
        end

        function computeStress(obj)
            strn  = permute(obj.strain.fValues,[1 3 2]);
            strn2(:,1,:,:) = strn;
            strs =squeeze(pagemtimes(obj.material.C,strn2));
            strs = permute(strs, [1 3 2]);

            z.mesh       = obj.mesh;
            z.fValues    = strs;
            z.quadrature = obj.quadrature;
            strFun = FGaussDiscontinuousFunction(z);

            obj.stress = strFun;
            obj.variables.stress = permute(strFun.fValues, [2 1 3]);
            obj.stressFun = strFun;
        end

        function computePrincipalDirection(obj)
            strss  = permute(obj.stressFun.fValues, [2 1 3]);
            s.type = obj.mesh.ndim;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(strss);
%             obj.variables.principalDirections = pcomp.direction;
%             obj.variables.principalStress     = pcomp.principalStress;
        end

    end

end
