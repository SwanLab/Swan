classdef MultigridTesting4 < handle

    properties (Access = public)
        u
    end

    properties (Access = private)
        nLevel
        multiLevelMesh

        ndim
        nDimf
        boundaryConditionsFine
        fineMaterial
        quad

        functionType
        nbasis

        Lt
        FineK

        FineKred
        coarseMesh
        CoarseK
        CoarseKred
        coarseMaterial
        boundaryConditionsCoarse
        fineMeshCoord
        fineMeshConnec
        fineMesh
        FineDispFun
        CoarseDispFun
        FineFred
        CoarseFred
        data


        I
        material
        boundaryConditions
        dispFun
        Kred
        Fred
        fem
    end

    methods (Access = public)

        function obj = MultigridTesting4()
            % tic
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init();
            obj.createMultiLevelMesh();
            mF           = obj.multiLevelMesh.mesh;
            coarseMeshes = obj.multiLevelMesh.coarseMeshes;
            interpolator = obj.multiLevelMesh.interpolator;

            s.type                = 'ELASTIC';
            s.scale               = 'MACRO';
            s.dim                 = '2D'; % Canviar per fer 2D
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'MULTIGRID';
            s.solverCase          = 'REDUCED';

            s.mesh                = mF;
            s.coarseMeshes        = coarseMeshes;
            s.interpolator        = interpolator;

            s.bc                  = obj.createBoundaryConditions(mF);
            ss.mesh               = mF;
            ss.boundaryConditions = s.bc;
            s.bcApplier           = BCApplier(ss);

            s.material            = obj.createMaterial(mF);
            s.dispFun             = LagrangianFunction.create(mF, obj.nDimf,obj.functionType);
            LHS                   = obj.computeStiffnessMatrix(mF,s.material,s.dispFun);
            RHS                   = obj.computeForces(mF,s.dispFun,s.bc,s.material,LHS,s.solverCase);
            s.LHS                 = s.bcApplier.fullToReducedMatrixDirichlet(LHS);
            s.RHS                 = s.bcApplier.fullToReducedVectorDirichlet(RHS);
            
            
            s.tol                 = 1e-6;
            s.nLevel              = obj.nLevel;
            s.nDimf               = obj.nDimf;
            tic
            solver                = Solver.create(s);
            obj.u                 = solver.solve();
            toc

            obj.postProcess();
            % obj.plotRes(obj.u,mF,s.bcApplier);
        end
    end

    methods (Access = private)

        function init(obj)
            obj.nDimf        = 2; % Canviar per fer 2D
            obj.nbasis       = 20;
            obj.functionType = 'P1';
            obj.nLevel       = 3;
            obj.ndim         = 2; % Canviar per fer 2D
        end

        function createMultiLevelMesh(obj)
            s.nX               = 2; % Canviar per fer 2D 
            s.nY               = 2; % Canviar per fer 2D 
            s.nZ               = 2; % Canviar per fer 2D 
            s.nLevel           = obj.nLevel;
            s.length           = 1; % Canviar per fer 2D
            s.height           = 1; % Canviar per fer 2D
            s.width            = 1; % Canviar per fer 2D
            s.ndim             = obj.ndim;
            m                  = MultilevelMesh(s);
            obj.multiLevelMesh = m;
        end

        function bc = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichlet = DirichletCondition(mesh,Dir);
            pointload = PointLoad(mesh,PL);
             % need this because force applied in the face not in a point
            pointload.values=pointload.values/size(pointload.dofs,1);
            fvalues = zeros(mesh.nnodes*obj.nDimf,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues = reshape(fvalues,obj.nDimf,[])';
            pointload.fun.setFValues(fvalues)

            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun =[];
            s.mesh = mesh;
            bc          = BoundaryConditions(s);
        end

        function dim = getFunDims(obj,mesh)
            s.fValues   = mesh.coord;
            s.mesh      = mesh;
            disp        = P1Function(s);
            d.ndimf     = disp.ndimf;
            d.nnodes    = size(disp.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = mesh.nnodeElem;
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir.domain    = @(coor) isLeft(coor);
            Dir.direction = [1,2,3];
            Dir.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -1;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E1  = 1;
            nu1 = 1/3;
            E   = ConstantFunction.create(E1,mesh);
            nu  = ConstantFunction.create(nu1,mesh);
            young   = E;
            poisson = nu;
        end

        function material = createMaterial(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function LHS = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            ndimf      = displacementFun.ndimf;
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = LagrangianFunction.create(mesh,ndimf, 'P1');
            s.trial    = displacementFun;
            s.material = material;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

         function forces = computeForces(obj,mesh,dispFun,boundaryConditions,material,stiffness,solverCase)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
%             s.dim       = obj.getFunDims();
            s.dim.ndofs = dispFun.nDofs;
            s.BC       = boundaryConditions;
            s.mesh     = mesh;
            s.material = material;
%             s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(solverCase,'REDUCED')
                R = RHSint.computeReactions(stiffness);
                forces = rhs+R;
            else
                forces = rhs;
            end
        end

        function postProcess(obj)
            DOFs = [1056 9312 25760 50400 83232];
            time2 = [0.051821 0.072445 0.113856 0.229354 0.282576];
            time5 = [0.00513 0.020559 0.046073 0.090326 0.130545];
            time10 = [0.004523 0.031546 0.067093 0.135769 0.16341];
            time20 = [0.005 0.046528 0.084207 0.174044 0.292117];

            figure(1)
            plot(DOFs,time2)
            hold on
            plot(DOFs,time5)
            plot(DOFs,time10)
            plot(DOFs,time20)
            title('Time vs DOFs')
            xlabel('DOFs')
            ylabel('Time(s)')

            legend('2 iters per level','5 iters per level','10 iters per level','20 iters per level','location','northwest')

            time2Total = [0.887724, 1.18978, 7.160798, 38.034618, 155.408853];
            time5Total = [0.33394, 1.047098, 7.062304, 39.194001, 163.584521];
            time10Total = [0.259715, 1.041644, 7.06799, 39.915648, 159.657016];
            time20Total = [0.266481, 0.97331, 7.010966, 38.85078, 151.357202];

            figure(2)
            plot(DOFs,time2Total)
            hold on
            plot(DOFs,time5Total)
            plot(DOFs,time10Total)
            plot(DOFs,time20Total)
            title('Total Time vs DOFs')
            xlabel('DOFs')
            ylabel('Time(s)')

            legend('2 iters per level','5 iters per level','10 iters per level','20 iters per level','location','northwest')

            DOFs3D = [1944, 13872, 45000, 104544, 201720, 345744, 545832, 753248];
            multigridTime = [0.863221, 14.131017, 159.585559, 916.336431, 3782.675412, 9980.822431, 24716.03475, 61790.087];
            directTime = [0.030516, 0.112095, 4.781604, 5.6041, 246.187, 3577.882947, 3760.532118, 109925.6268];

            figure(3)
            plot(DOFs3D,multigridTime)
            hold on
            plot(DOFs3D,directTime)
            title('Total Time vs DOFs')
            xlabel('DOFs')
            ylabel('Time(s)')

            legend('Multigrid', 'Direct','Location','northwest')
        end

        function plotRes(obj,res,mesh,bcApplier,numItr)
            xFull = bcApplier.reducedToFullVectorDirichlet(res);
            s.fValues = reshape(xFull,obj.nDimf,[])';
            s.mesh = mesh;
            if obj.nDimf<3
                s.fValues(:,end+1) = 0;
            end
            s.order   = 'P1';
%             s.ndimf = obj.nDimf;
            uFeFun = LagrangianFunction(s);
%             xF = P1Function(s);
            %xF.plot();
            uFeFun.print('uPrueva','Paraview')
            uFeFun.print('uPrueva','Paraview')
            fclose('all');
        end
    end

end