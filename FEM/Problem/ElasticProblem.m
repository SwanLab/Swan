classdef ElasticProblem < handle
    
    properties (Access = public)
        variables
    end

    properties (Access = private)
        boundaryConditions
        displacement
        problemData
        stiffnessMatrix
        RHS
        solver
        geometry
    end

    properties (Access = protected)
        quadrature
        dim, dim2
        material

        vstrain

        mesh % For Homogenization
    end

    methods (Access = public)

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.computeDimensions();
            obj.computeNewDimensions();
            obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
%             obj.computeStiffnessMatrixOld();
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStrain();
            obj.computeStress();
            obj.computePrincipalDirection();
        end

        function plot(obj)
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.dim;
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
       
        function print(obj,filename)
            s.quad = obj.quadrature;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.createVariablesToPrint();
            s.ptype     = obj.problemData.ptype;
            s.ndim      = obj.dim.ndim;
            s.pdim      = obj.problemData.pdim;
            s.type      = obj.createPrintType();
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            pd.scale        = cParams.scale;
            pd.pdim         = cParams.dim;
            pd.ptype        = cParams.type;
            pd.bc.dirichlet = cParams.bc.dirichlet;
            pd.bc.pointload = cParams.bc.pointload;
            obj.problemData = pd;
            obj.createQuadrature();
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function computeDimensions(obj)
            s.ngaus = obj.quadrature.ngaus;
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function computeNewDimensions(obj)
            s.fieldName = 'u';
            s.mesh = obj.mesh;
            s.ndimf = str2double(regexp(obj.problemData.pdim,'\d*','Match'));

            d = DimensionVector(s);
            d.create(s)
            d.applyNgaus(obj.quadrature.ngaus);
            obj.dim = d;
        end

        function createBoundaryConditions(obj)
            s.dim        = obj.dim;
            s.mesh       = obj.mesh;
            s.scale      = obj.problemData.scale;
            s.bc         = obj.problemData.bc;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrix';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            s.material     = obj.material;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            obj.stiffnessMatrix = K;
        end

        function computeStiffnessMatrixOld(obj)
            s.type = 'ElasticStiffnessMatrixOld';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            s.material     = obj.material;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            obj.stiffnessMatrix = K;
        end

        function computeForces(obj)
            s.dim         = obj.dim;
            s.BC          = obj.boundaryConditions;
            s.mesh        = obj.mesh;
            s.material    = obj.material;
            s.globalConnec = obj.mesh.connec;
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            fcomp = ForcesComputer(s);
            f = fcomp.compute();
            R = fcomp.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = f + R;
            obj.RHS = f;
        end

        function u = computeDisplacements(obj)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.stiffnessMatrix);
            Fred = bc.fullToReducedVector(obj.RHS);
            u = obj.solver.solve(Kred,Fred);
            u = bc.reducedToFullVector(u);
            obj.variables.d_u = u;
        end

        function computeStrain(obj)
            s.dim                = obj.dim;
            s.mesh               = obj.mesh;
            s.displacement       = obj.variables.d_u;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;
        end

        function computeStress(obj)
            s.C      = obj.material.C;
            s.dim    = obj.dim;
            s.strain = obj.variables.strain;
            scomp  = StressComputer(s);
            stress = scomp.compute();
            obj.variables.stress = stress;
        end

        function computePrincipalDirection(obj)
            stress = obj.variables.stress;
            s.type = obj.problemData.pdim;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(stress);
            obj.variables.principalDirections = pcomp.direction;
            obj.variables.principalStress     = pcomp.principalStress;
        end

        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh    = obj.mesh;
            dI.outName = fileName;
            dI.pdim = '2D';
            dI.ptype = 'ELASTIC';
            ps = PostProcessDataBaseCreator(dI);
            d = ps.getValue();
        end

        function uM = splitDisplacement(obj)
            u = obj.variables.d_u;
            nu = obj.dim.ndimField;
            nnode = round(length(u)/nu);
            nodes = 1:nnode;
            uM = zeros(nnode,nu);
            for idim = 1:nu
                dofs = nu*(nodes-1)+idim;
                uM(:,idim) = u(dofs);
            end
        end


    end

    methods (Access = protected)

        function f = createVariablesToPrint(obj)
            f = obj.variables;
            f.u = obj.splitDisplacement();
        end

        function t = createPrintType(obj)
           t = 'Elasticity';
        end


    end

end
