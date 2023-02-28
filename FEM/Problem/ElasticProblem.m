classdef ElasticProblem < handle
    
    properties (Access = public)
        variables
        boundaryConditions
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        LHS
        RHS
        solver
        geometry
        scale
        pdim
        ptype
        inputBC

        strain
        stress
    end

    properties (Access = protected)
        quadrature
        material

        vstrain

        mesh % For Homogenization
        interpolation
        interpTranslator
        interpolationType
        displacementField
        displacementFun
    end

    methods (Access = public)

        function obj = ElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementField();
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
            s.displacement = obj.variables.d_u;
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
       
        function print(obj,filename)
            [fun, funNames] = obj.getFunsToPlot();
            a.mesh     = obj.mesh;
            a.filename = filename;
            a.fun      = fun;
            a.funNames = funNames;
            pst = ParaviewPostprocessor(a);
            pst.print();
        end

        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.uFun{:}, obj.strainFun{:}, obj.stressFun{:}};
            funNames = {'displacement', 'strain', 'stress'};
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
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end

        function createDisplacementField(obj)
            ndimf = regexp(obj.pdim,'\d*','Match');
            s.mesh               = obj.mesh;
            s.ndimf              = str2double(ndimf);
            s.interpolationOrder = obj.interpolationType; %obj.interpolationType
            f = Field(s);
            obj.inputBC = f.translateBoundaryConditions(obj.inputBC);
            obj.displacementField = f;


            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            nNodes = size(obj.mesh.coord,1);
            s.ndimf   = nDimf;
            s.mesh    = obj.mesh;
            s.fValues = zeros(nNodes,nDimf);
            f = P1Function(s);
            obj.displacementFun = f;
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
            s.type =  'DIRECT';
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
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS);
            obj.variables.fext = rhs + R;
            obj.RHS = rhs;
        end

        function u = computeDisplacements(obj)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.LHS);
            Fred = bc.fullToReducedVector(obj.RHS);
            u = obj.solver.solve(Kred,Fred);
            u = bc.reducedToFullVector(u);
            obj.variables.d_u = u;

            z.mesh   = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFeFun = P1Function(z);
            obj.uFun{end+1} = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrain(obj)
            s.dim.ndimf    = obj.displacementFun.ndimf;
            s.mesh         = obj.mesh;
            s.quadrature   = obj.quadrature;
            s.displacement = obj.variables.d_u;
            s.dispField    = obj.displacementField;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;
            
%             strFun = obj.uFun.computeSymmetricGradient(obj.quadrature);
%             strFun.applyVoigtNotation();
            z.mesh       = obj.mesh;
            z.fValues    = permute(strain, [2 1 3]);
            z.quadrature = obj.quadrature;
            strFun = FGaussDiscontinuousFunction(z);

            wtfun = obj.displacementFun.computeSymmetricGradient(obj.quadrature);
            wtfun.applyVoigtNotation();

            perm = permute(wtfun.fValues, [2 1 3]);
            norm(squeeze(perm(:)-strain(:)))/norm(strain(:))
            obj.variables.strain = perm;
            obj.strainFun{end+1} = wtfun;
%             strnFun = obj.displacementFun.computeSymmetricGradient(obj.quadrature);
%             strnFun.applyVoigtNotation();
%             obj.strain = strnFun;
%             obj.strainFun{end+1} = strnFun;
%             obj.variables.strain = permute(strnFun.fValues, [2 1 3]);
        end

        function computeStress(obj)
            s.C      = obj.material.C;
            s.strain = obj.variables.strain;
            scomp  = StressComputer(s);
            stress = scomp.compute();
            obj.variables.stress = stress;
            
            z.mesh       = obj.mesh;
            z.fValues    = permute(stress, [2 1 3]);
            z.quadrature = obj.quadrature;
            strFun = FGaussDiscontinuousFunction(z);

            obj.stressFun{end+1} = strFun;

%             strn  = permute(obj.strain.fValues,[1 3 2]);
%             strn2(:,1,:,:) = strn;
%             strs =squeeze(pagemtimes(obj.material.C,strn2));
%             strs = permute(strs, [1 3 2]);
%             z.mesh       = obj.mesh;
%             z.fValues    = strs;
%             z.quadrature = obj.quadrature;
%             stressF = FGaussDiscontinuousFunction(z);
%             obj.stress = stressF;
%             
%             obj.variables.stress = permute(stressF.fValues, [2 1 3]);
%             obj.stressFun{end+1} = stressF;
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

    end

end
