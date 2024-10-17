classdef PreconditionerDirichletNeumann < handle
    
    properties (Access = public)

    end

    properties (Access = private)
        subdomainLHS
        isDirichlet
    end

    properties (Access = private)
        LHS
        bcApplier
        DirCond
        ddDofManager
        RHS       
        weight
        subdomainMesh
        LHSreduced
    end

    methods (Access = public)

        function obj = PreconditionerDirichletNeumann(cParams)
            obj.init(cParams);
            obj.createSubdomainsBCappliers();            
            obj.computeSubdomainLHS();
            obj.computeSubdomainLHSreduced();
         %   obj.DirichletNeumannSolver();

            %

            % %             obj.obtainCornerNodes();
            % %             fineMesh = MeshFromRVE
            %             obj.createSubDomainMeshes();
            %             obj.createInterfaceSubDomainMeshes();
            %             obj.createDomainMesh();


            %             s.referenceMesh = obj.referenceMesh;
            %             mC = MeshCreatorFromSubmeshes();
            %             obj.meshDomain = mC.mesh;

            %             preconditioner = obj.createPreconditioner(mC.submeshes);

        end

        function z = apply(obj,r)

            z = 0.3*r;
            lhsS = obj.LHS;
            rR = obj.bcApplier.reducedToFullVectorDirichlet(r);
            rS = obj.ddDofManager.global2local(rR);

        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ddDofManager  = cParams.ddDofManager;
            obj.bcApplier     = cParams.bcApplier;
            obj.weight        = 0.5;
            obj.DirCond       = cParams.DirCond;
            obj.LHS           = cParams.LHS;
            obj.subdomainMesh = cParams.subdomainMesh;
            obj.isDirichlet{1} = false;
            obj.isDirichlet{2} = true;
        end

        function createSubdomainsBCappliers(obj)
            for iS = 1:numel(obj.isDirichlet)
                 s.mesh                  = obj.subdomainMesh{iS};
                 s.boundaryConditions    = obj.boundaryConditions;
                 obj.bcApplier           = BCApplier(s);                

            end
        end


        function computeSubdomainLHS(obj)
            sLHS             = obj.ddDofManager.global2localMatrix(obj.LHS);
            obj.subdomainLHS = obj.ddDofManager.scaleInterfaceValuesMatrix(sLHS,obj.weight);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isBottom  = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop  = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor)| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = 2;
            PL.value     = -0.1;
        end

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.fValues   = fvalues;

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end


        function DirichletNeumannSolver(obj)
            tol=1e-8;
            e=1;
            theta = 3/4;
            iter=1;
            while e>tol
                %                 for i=1:obj.nSubdomains(1)
                obj.computeForces(2);
                Kred = obj.boundaryConditions{2}.fullToReducedMatrix(obj.LHS{2});
                Fred = obj.boundaryConditions{2}.fullToReducedVector(obj.RHS{2});
                u{2} = Kred\Fred;
                u{2} = obj.boundaryConditions{2}.reducedToFullVector(u{2});
                R = obj.RHS{2}-0*obj.LHS{2}*u{2};
                obj.plotSolution(R,obj.meshSubDomain{2}, 2,iter)
                %                 obj.plotSolution(u{2},obj.meshSubDomain{2}, 2,iter)

                obj.computeForces(1);
                R_int = obj.LHS{2}(obj.interfaceDof(:,2),:)*u{2};
                obj.RHS{1}(obj.interfaceDof(:,1)) = obj.RHS{1}(obj.interfaceDof(:,1)) - R_int;

                Kred = obj.boundaryConditions{1}.fullToReducedMatrix(obj.LHS{1});
                Fred = obj.boundaryConditions{1}.fullToReducedVector(obj.RHS{1});
                u{1} = Kred\Fred;
                u{1} = obj.boundaryConditions{1}.reducedToFullVector(u{1});
                R = obj.RHS{1}-0*obj.LHS{1}*u{1};
                obj.plotSolution(R,obj.meshSubDomain{1}, 1,iter)
                %                 obj.plotSolution(u{1},obj.meshSubDomain{1}, 1,iter)

                u1int = theta*u{2}(obj.interfaceDof(:,2)) + (1-theta)*u{1}(obj.interfaceDof(:,1));

                for idof = 1: length(obj.interfaceDof(:,2))
                    ind = find(obj.boundaryConditions{2}.dirichlet == obj.interfaceDof(idof,2));
                    obj.boundaryConditions{2}.dirichlet_values(ind) = u1int(idof);
                end

                e = norm(u1int-u{2}(obj.interfaceDof(:,2)));

                iter=iter+1;
                %                 end
            end
        end

        function plotSolution(obj,x,mesh,domain,iter)
            %             xFull = bc.reducedToFullVector(x);
            s.fValues = reshape(x,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %     xF.plot();
            xF.print(['domain',num2str(domain),'_',num2str(iter)],'Paraview')
            fclose('all');
        end

    end
    
end