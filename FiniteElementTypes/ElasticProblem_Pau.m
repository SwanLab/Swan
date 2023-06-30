classdef ElasticProblem_Pau < handle
    
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
        scale
        pdim
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

        function obj = ElasticProblem_Pau(cParams)
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
            obj.inputBC     = cParams.bc;
            if isprop(cParams, 'interpolationType') % later on for P2
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end
            obj.createQuadrature();
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('ORDER10');
            obj.quadrature = quad;
        end

        function createDisplacementFun(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
%             obj.displacementFun = P1Function.create(obj.mesh, nDimf);
            k=2;
            l=LagrangeElement.create("SIMPLICIAL",k,2);
            obj.displacementFun = FE_LagrangianFunction.create(obj.mesh,nDimf,k,l);
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
%             d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.nnodeElem = size(obj.displacementFun.interpolation.finiteElement.shapeFunctions,1);
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBoundaryConditions(obj)
            dim = obj.getFunDims();
            bc = obj.inputBC;
bc.dirichlet = [
    1, 1, 0;
    1, 2, 0;
    3, 1, 0;
    3, 2, 0;
    8, 1, 0;
    8, 2, 0;
    15, 1, 0;
    15, 2, 0;
    22, 1, 0;
    22, 2, 0;
    43, 1, 0;
    43, 2, 0;
    49, 1, 0;
    49, 2, 0;
    53, 1, 0;
    53, 2, 0;
    75, 1, 0;
    75, 2, 0;
    147, 1, 0;
    147, 2, 0;
    153, 1, 0;
    153, 2, 0;
    157, 1, 0;
    157, 2, 0;
    179, 1, 0;
    179, 2, 0;
    208, 1, 0;
    208, 2, 0;
    209, 1, 0;
    209, 2, 0;
    246, 1, 0;
    246, 2, 0;
    247, 1, 0;
    247, 2, 0;
    547, 1, 0;
    547, 2, 0;
    553, 1, 0;
    553, 2, 0;
    557, 1, 0;
    557, 2, 0;
    579, 1, 0;
    579, 2, 0;
    608, 1, 0;
    608, 2, 0;
    609, 1, 0;
    609, 2, 0;
    646, 1, 0;
    646, 2, 0;
    647, 1, 0;
    647, 2, 0;
    760, 1, 0;
    760, 2, 0;
    761, 1, 0;
    761, 2, 0;
    792, 1, 0;
    792, 2, 0;
    793, 1, 0;
    793, 2, 0;
    814, 1, 0;
    814, 2, 0;
    815, 1, 0;
    815, 2, 0;
    940, 1, 0;
    940, 2, 0;
    941, 1, 0;
    941, 2, 0;
    2115, 1, 0;
    2115, 2, 0;
    2121, 1, 0;
    2121, 2, 0;
    2125, 1, 0;
    2125, 2, 0;
    2147, 1, 0;
    2147, 2, 0;
    2176, 1, 0;
    2176, 2, 0;
    2177, 1, 0;
    2177, 2, 0;
    2214, 1, 0;
    2214, 2, 0;
    2215, 1, 0;
    2215, 2, 0;
    2328, 1, 0;
    2328, 2, 0;
    2329, 1, 0;
    2329, 2, 0;
    2360, 1, 0;
    2360, 2, 0;
    2361, 1, 0;
    2361, 2, 0;
    2382, 1, 0;
    2382, 2, 0;
    2383, 1, 0;
    2383, 2, 0;
    2508, 1, 0;
    2508, 2, 0;
    2509, 1, 0;
    2509, 2, 0;
    2920, 1, 0;
    2920, 2, 0;
    2921, 1, 0;
    2921, 2, 0;
    2952, 1, 0;
    2952, 2, 0;
    2953, 1, 0;
    2953, 2, 0;
    2974, 1, 0;
    2974, 2, 0;
    2975, 1, 0;
    2975, 2, 0;
    3100, 1, 0;
    3100, 2, 0;
    3101, 1, 0;
    3101, 2, 0;
    3254, 1, 0;
    3254, 2, 0;
    3255, 1, 0;
    3255, 2, 0;
    3258, 1, 0;
    3258, 2, 0;
    3259, 1, 0;
    3259, 2, 0;
    3474, 1, 0;
    3474, 2, 0;
    3475, 1, 0;
    3475, 2, 0;
    3478, 1, 0;
    3478, 2, 0;
    3479, 1, 0;
    3479, 2, 0
];
            
% bc.dirichlet = [
%     1 1 0;
%     1 2 0;
%     3 1 0;
%     3 2 0;
%     8 1 0;
%     8 2 0;
%     15 1 0;
%     15 2 0;
%     22 1 0;
%     22 2 0;
%     43 1 0;
%     43 2 0;
%     49 1 0;
%     49 2 0;
%     53 1 0;
%     53 2 0;
%     75 1 0;
%     75 2 0;
%     147 1 0;
%     147 2 0;
%     153 1 0;
%     153 2 0;
%     157 1 0;
%     157 2 0;
%     179 1 0;
%     179 2 0;
%     208 1 0;
%     208 2 0;
%     209 1 0;
%     209 2 0;
%     246 1 0;
%     246 2 0;
%     247 1 0;
%     247 2 0
% ];
            bc.ndimf = dim.ndimf;
            bc.ndofs = dim.ndofs;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
%             bc.dirichlet = [1;2;5;6;15;16;29;30;43;44];
%             bc.dirichlet_values = [0;0;0;0;0;0;0;0;0;0];
%             bc.dirichlet = [1;2;5;6;15;16;31;32;33;34;59;60;61;62];
%             bc.dirichlet_values = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
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
            s.quadratureOrder = 'ORDER10';
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

            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.displacementFun.ndofs])';
            z.polynomialOrder = 2;
            z.finiteElement = LagrangeElement.create("SIMPLICIAL",z.polynomialOrder,2);
            uFeFun = FE_LagrangianFunction(z);
            obj.uFun{end+1} = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.displacementFun.ndofs])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrain(obj)
            strFun = obj.displacementFun.computeSymmetricGradient(obj.quadrature);
            strFun.applyVoigtNotation();
            perm = permute(strFun.fValues, [2 1 3]);
            obj.variables.strain = perm;
            obj.strainFun{end+1} = strFun;
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
            obj.stressFun{end+1} = strFun;
        end

        function computePrincipalDirection(obj)
            strss  = obj.variables.stress;
            s.type = obj.pdim;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(strss);
            obj.variables.principalDirections = pcomp.direction;
            obj.variables.principalStress     = pcomp.principalStress;
        end

    end

end
