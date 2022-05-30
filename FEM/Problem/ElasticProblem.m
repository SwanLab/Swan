classdef ElasticProblem < handle
    
    properties (Access = public)
        variables
    end

    properties (Access = private)
        boundaryConditions
        displacement
        stiffnessMatrix
        RHS
        solver
        geometry

        scale
        pdim
        ptype
        inputBC
    end

    properties (Access = protected)
        quadrature
        dim
        material

        vstrain

        mesh % For Homogenization
        interpolation
        interpTranslator
        interpolationType
        displacementField
    end

    methods (Access = public)

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.computeDimensions();
            obj.createDisplacementField();
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
            s.dim          = obj.displacementField.dim;
            s.mesh         = obj.mesh;
            s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.displacementField.dim;
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
            s.ptype     = obj.ptype;
            s.ndim      = obj.displacementField.dim.ndimf;
            s.pdim      = obj.pdim;
            s.type      = obj.createPrintType();
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.scale       = cParams.scale;
            obj.pdim        = cParams.dim;
            obj.ptype       = cParams.type;
            obj.inputBC     = cParams.bc;
            if isprop(cParams, 'interpolationType')
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end
            obj.createQuadrature();
            obj.createInterpolation();
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,obj.interpolationType);
%             int = Interpolation.create(obj.mesh,'QUADRATIC');
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end

        function computeDimensions(obj)
            s.type      = 'Vector';
            s.fieldName = 'u';
            s.mesh      = obj.mesh;
            s.ndimf = str2double(regexp(obj.pdim,'\d*','Match'));
            d = DimensionVariables.create(s);
            d.compute()
            obj.dim = d;
        end

        function createDisplacementField(obj)
            ndimf = regexp(obj.pdim,'\d*','Match');
            s.mesh               = obj.mesh;
            s.ndimf              = str2double(ndimf);
            s.inputBC            = obj.inputBC;
            s.scale              = obj.scale;
            s.interpolationOrder = obj.interpolationType; %obj.interpolationType
            obj.displacementField = Field(s);
        end

        function createBoundaryConditions(obj)
            bc = obj.displacementField.inputBC;
            bc.ndimf = obj.displacementField.dim.ndimf;
            bc.ndofs = obj.displacementField.dim.ndofs;
            s.dim   = obj.displacementField.dim;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = obj.displacementField.dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrix';
            s.mesh          = obj.mesh;
            s.globalConnec  = obj.displacementField.connec;
            s.dim           = obj.displacementField.dim;
            s.material      = obj.material;
            s.interpolation = obj.interpolation;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            obj.stiffnessMatrix = K;
        end

        function computeStiffnessMatrixOld(obj)
            s.type = 'ElasticStiffnessMatrixOld';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.displacementField.connec;
            s.dim          = obj.displacementField.dim;
            s.material     = obj.material;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            obj.stiffnessMatrix = K;
        end

        function computeForces(obj)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.displacementField.dim;
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.displacementField.connec;
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = rhs + R;
            obj.RHS = rhs;
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
            s.dim          = obj.displacementField.dim;
            s.mesh         = obj.mesh;
            s.quadrature   = obj.quadrature;
            s.displacement = obj.variables.d_u;
            s.dispField    = obj.displacementField;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;
        end

        function computeStress(obj)
            s.C      = obj.material.C;
            s.dim    = obj.displacementField.dim;
            s.strain = obj.variables.strain;
            scomp  = StressComputer(s);
            stress = scomp.compute();
            obj.variables.stress = stress;
        end

        function computePrincipalDirection(obj)
            stress = obj.variables.stress;
            s.type = obj.pdim;
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
            nu = obj.displacementField.ndimf;
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
