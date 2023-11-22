classdef ElasticProblem < handle
    
    properties (Access = public)
%         variables
        boundaryConditions
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        LHS
        RHS
        solver
        scale
        inputBC
        strain
        stress
    end

    properties (Access = protected)
        quadrature
        material
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
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createSolver(obj)
            s.type =  'rMINRES';
            obj.solver = Solver.create(s);
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacementFun;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeForces(obj)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
%             s.globalConnec = obj.displacementField.connec;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS);
%             obj.variables.fext = rhs + R;
            obj.RHS = rhs;
        end

        function u = computeDisplacements(obj)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.LHS);
            Fred = bc.fullToReducedVector(obj.RHS);
            u = obj.solver.solve(Kred,Fred);
            u = bc.reducedToFullVector(u);
%             obj.variables.d_u = u;

            z.mesh    = obj.mesh;
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
%             obj.variables.strain = perm;
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
%             obj.variables.stress = permute(strFun.fValues, [2 1 3]);
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
