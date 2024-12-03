classdef DirichletNeumann < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        meshDomain
        boundaryConditions
        material
        meshReference
        interfaceMeshReference
        meshSubDomain
        interfaceConnec
        nSubdomains
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        cornerNodes
        quad
        interfaceDof

        displacementFun
        LHS
        RHS
        scale
    end

    methods (Access = public)

        function obj = DirichletNeumann()
            close all
            obj.init();
            obj.createReferenceMesh();

            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec] = m.create();


            obj.createDisplacementFun();
            obj.interfaceDof = obj.computeLocalInterfaceDof;

            for idom=1:obj.nSubdomains(1)
                obj.quad{idom} = Quadrature.set(obj.meshSubDomain{idom}.type);
                obj.quad{idom}.computeQuadrature('QUADRATIC');
                obj.createDomainMaterial(idom);
                obj.computeStiffnessMatrix(idom);
                rawBC = obj.createRawBoundaryConditions(idom);
                obj.boundaryConditions{idom} = obj.createBoundaryConditions(obj.meshSubDomain{idom},rawBC{idom},idom);
            end

            obj.DirichletNeumannSolver();

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

            obj.solveDomainProblem();



        end

        function solveDomainProblem(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverTyp = 'PCG';
            s.preconditionerType = 'EIFEM';
            s.tol = 1e-6;

            fem        = FEM.create(s);
            fem.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [2 1]; %nx ny
            obj.scale    = 'MACRO';
        end

        function createReferenceMesh(obj)
            %             filename   = 'lattice_ex1';
            %             a.fileName = filename;
            %             femD       = FemDataContainer(a);
            %             mS         = femD.mesh;
            %             bS         = mS.createBoundaryMesh();
            % Generate coordinates
            x1 = linspace(0,1,10);
            x2 = linspace(0,1,10);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord    = coord(:,1:2);
            s.connec   = F;
            mS         = Mesh(s);
            bS         = mS.createBoundaryMesh();

            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function cMesh = createCoarseMesh(obj)
            xmax = max(obj.meshReference.coord(:,1));
            xmin = min(obj.meshReference.coord(:,1));
            ymax = max(obj.meshReference.coord(:,2));
            ymin = min(obj.meshReference.coord(:,2));
            coord(1,1) = xmax;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymax;
            coord(3,1) = xmin;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymin;
            connec = [1 2 3 4];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh(s);
        end

        function createDomainMaterial(obj,i)
            %             ngaus = 1;
            %             m = obj.meshSubDomain{i};
            obj.material{i} = obj.createMaterial(i);
        end

        function L = computeReferenceMeshLength(obj)
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1));
            Ly = max(coord(:,2));
            L = [Lx Ly];
        end

        function BC = createRawBoundaryConditions(obj,i)

            dirichletBc{1}.boundaryId=[1];
            dirichletBc{1}.dof{1}=[1,2];
            dirichletBc{1}.value{1}=[0,0];
            %             newmanBc{1}.boundaryId=[4,1];
            %             newmanBc{1}.dof{1}=2;
            %             newmanBc{1}.dof{2}=[1,2];
            %             newmanBc{1}.value{1}=10;
            %             newmanBc{1}.value{2}=[0,0];
            newmanBc{1}=[];

            dirichletBc{2}.boundaryId=[1];
            dirichletBc{2}.dof{1}=[1,2];
            %             dirichletBc{2}.dof{2}=[1,2];
            dirichletBc{2}.value{1}=[0,0];
            %             dirichletBc{2}.value{2}=[0,0];
            newmanBc{2}.boundaryId=2;
            newmanBc{2}.dof{1}=[2];
            newmanBc{2}.value{1}=[-1];



            nx= obj.nSubdomains(1);
            %             for i=1:nx
            bM = obj.meshSubDomain{i}.createBoundaryMesh();
            [dirichlet,pointload] = obj.createBc(bM,dirichletBc{i},newmanBc{i});
            BC{i}.dirichlet=dirichlet;
            BC{i}.pointload=pointload;
            %             end


            %             obj.boundaryConditions = BC;
        end

        function bc = createBoundaryConditions(obj,mesh,bcV,i)
            dim = obj.getFunDims(i);
            bcV.ndimf = dim.ndimf;
            bcV.ndofs = dim.ndofs;
            s.mesh  = mesh;
            s.scale = 'MACRO';
            s.bc    = {bcV};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
        end

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)

            dirichlet = obj.createBoundaryCondition(boundaryMesh,dirchletBc);
            if ~isempty(newmanBc)
                pointload = obj.createBoundaryCondition(boundaryMesh,newmanBc);
            else
                pointload=[];
            end
        end

        function cond = createBoundaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                %                 ncond  = length(condition.dof(nbound,:));
                ncond  = length(condition.dof{ibound});
                nodeId= reshape(unique(bM{condition.boundaryId(ibound)}.globalConnec),[],1);
                nbd   = length(nodeId);
                condition.value{ibound} = condition.value{ibound}/nbd;
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof{ibound}(icond),[nbd,1]), repmat(condition.value{ibound}(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end

        function material = createMaterial(obj,i)
            I = ones(obj.meshSubDomain{i}.nelem,obj.quad{i}.ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.meshSubDomain{i}.nelem;
            s.mesh  = obj.meshSubDomain{i};
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            mat = Material.create(s);
            mat.compute(s);
            material = mat;
        end

        function createDisplacementFun(obj)
            nx= obj.nSubdomains(1);
            for i=1:nx
                obj.displacementFun{i} = P1Function.create(obj.meshSubDomain{i}, obj.meshSubDomain{i}.ndim);
            end
        end

        function dim = getFunDims(obj,i)
            d.ndimf  = obj.displacementFun{i}.ndimf;
            d.nnodes = size(obj.displacementFun{i}.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function computeStiffnessMatrix(obj,i)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshSubDomain{i};
            s.fun      = obj.displacementFun{i};
            s.material = obj.material{i};
            lhs = LHSintegrator.create(s);
            obj.LHS{i} = lhs.compute();
        end

        function computeForces(obj,i)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims(i);
            s.BC       = obj.boundaryConditions{i};
            s.mesh     = obj.meshSubDomain{i};
            s.material = obj.material{i};
            %             s.globalConnec = obj.displacementFun{i}.connec;
            s.globalConnec = obj.meshSubDomain{i}.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS{i});
            %             obj.variables.fext = rhs + R;
            obj.RHS{i} = rhs+R;
        end

        function interfaceDof = computeLocalInterfaceDof(obj)
            nint = size(obj.interfaceConnec,2);
            globaldof=0;
            for i=1:nint
                dofaux=0;
                dim = obj.getFunDims(i);
                nodesI = obj.interfaceConnec(:,i);

                for iunkn=1:dim.ndimf
                    DOF = dim.ndimf*(nodesI - 1) + iunkn;
                    DOF = DOF-globaldof;
                    dofaux= [dofaux; DOF];
                end
                interfaceDof(:,i) = dofaux(2:end);
                globaldof = globaldof + dim.ndofs;
            end
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
