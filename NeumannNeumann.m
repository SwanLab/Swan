classdef NeumannNeumann < handle

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
        theta
    end

    methods (Access = public)

        function obj = NeumannNeumann()
            close all
            obj.init();
            obj.createReferenceMesh();

            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec] = m.create();

            
            obj.createDisplacementFun();
            obj.interfaceDof = obj.computeLocalInterfaceDof;
            obj.computeSubdomainLHS();
            
            obj.createSubdomainBoundaryConditions();

            obj.NeumannNeumannSolver();

         
%             obj.solveDomainProblem();



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
            obj.theta=0.5;
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
      
        function createDomainMaterial(obj,i)
            %             ngaus = 1;
            %             m = obj.meshSubDomain{i};
            obj.material{i} = obj.createMaterial(i);
        end
     
        function BC = createRawBoundaryConditionsDirichlet(obj,i)

            dirichletBc{1}.boundaryId=[1,2];
            dirichletBc{1}.dof{1}=[1,2];
            dirichletBc{1}.dof{2}=[1,2];
            dirichletBc{1}.value{1}=[0,0];
            dirichletBc{1}.value{2}=[0,0];

%             newmanBc{1}=[];
            newmanBc{1}.boundaryId=4;
            newmanBc{1}.dof{1}=[2];
            newmanBc{1}.value{1}=[-1];

            dirichletBc{2}.boundaryId=[1,2];
            dirichletBc{2}.dof{1}=[1,2];
            dirichletBc{2}.dof{2}=[1,2];
            dirichletBc{2}.value{1}=[0,0];
            dirichletBc{2}.value{2}=[0,0];

            newmanBc{2}.boundaryId=4;
            newmanBc{2}.dof{1}=[2];
            newmanBc{2}.value{1}=[-1];

            bM = obj.meshSubDomain{i}.createBoundaryMesh();
            [dirichlet,pointload] = obj.createBc(bM,dirichletBc{i},newmanBc{i});
            BC.dirichlet=dirichlet;
            BC.pointload=pointload;

        end

        function BC = createRawBoundaryConditionsNeumann(obj,i)

            dirichletBc{1}.boundaryId=[1];
            dirichletBc{1}.dof{1}=[1,2];
            dirichletBc{1}.value{1}=[0,0];

%             newmanBc{1}=[];
            newmanBc{1}.boundaryId=[2];
%             newmanBc{1}.dof{1}=[2];
%             newmanBc{1}.value{1}=[-1];
            newmanBc{1}.dof{1}=[1,2];
            newmanBc{1}.value{1}=[0,0];

            dirichletBc{2}.boundaryId=[2];
            dirichletBc{2}.dof{1}=[1,2];
            dirichletBc{2}.value{1}=[0,0];

            newmanBc{2}.boundaryId=[1];
%             newmanBc{2}.dof{1}=[2];
%             newmanBc{2}.value{1}=[-1];
            newmanBc{2}.dof{1}=[1,2];
            newmanBc{2}.value{1}=[0,0];


            bM = obj.meshSubDomain{i}.createBoundaryMesh();
            [dirichlet,pointload] = obj.createBc(bM,dirichletBc{i},newmanBc{i});
            BC.dirichlet=dirichlet;
            BC.pointload=pointload;

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

        function createSubdomainBoundaryConditions(obj)
            for idom=1:obj.nSubdomains(1)
                rawBC = obj.createRawBoundaryConditionsDirichlet(idom);
                obj.boundaryConditions.dirichletStep{idom} = obj.createBoundaryConditions(obj.meshSubDomain{idom},rawBC,idom);
            end
            
            for idom=1:obj.nSubdomains(1)
                rawBC = obj.createRawBoundaryConditionsNeumann(idom);
                obj.boundaryConditions.neumannStep{idom} = obj.createBoundaryConditions(obj.meshSubDomain{idom},rawBC,idom);
            end
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
            nx = obj.nSubdomains(1);
            ny = obj.nSubdomains(2);
            for i=1:ny
                for j=1:nx
                    obj.displacementFun{i,j} = P1Function.create(obj.meshSubDomain{i,j}, obj.meshSubDomain{i,j}.ndim);
                end
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

        function computeForces(obj,i,boundaryConditions)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims(i);
            s.BC       = boundaryConditions{i};
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

        function computeSubdomainLHS(obj)
            for idom=1:obj.nSubdomains(1)
                obj.quad{idom} = Quadrature.set(obj.meshSubDomain{idom}.type);
                obj.quad{idom}.computeQuadrature('QUADRATIC');
                obj.createDomainMaterial(idom);
                obj.computeStiffnessMatrix(idom);
            end
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

        function NeumannNeumannSolver(obj)
            tol=1e-8;
            e=1;
            theta = 3/4;
            iter=1;
            uInt = zeros(size(obj.interfaceDof,1),1); 
            while e>tol

                uD = obj.dirichletStep();
                obj.updateNeumanStepBC(uD);
                uN = obj.neumannStep();
                uInt_new = obj.updateDirichletStepBC(uN);

                e = norm(uInt_new-uInt)/norm(uInt);
                uInt = uInt_new;

                iter=iter+1;
                %                 end
            end
        end

        function u = dirichletStep(obj)
            for idom = 1:obj.nSubdomains(1)
                obj.computeForces(idom,obj.boundaryConditions.dirichletStep)
                Kred    = obj.boundaryConditions.dirichletStep{idom}.fullToReducedMatrix(obj.LHS{idom});
                Fred    = obj.boundaryConditions.dirichletStep{idom}.fullToReducedVector(obj.RHS{idom});
                u{idom} = Kred\Fred;
            end
        end

        function R = updateNeumanStepBC(obj,u)
            R = obj.computeInterfaceReactions(u);
            obj.updateNeumanValues(R);
        end

        function R = computeInterfaceReactions(obj,u)   
             w = [obj.theta,1- obj.theta];
             R=zeros(size(obj.interfaceDof,1),1);
             for idom=1:obj.nSubdomains(1)
                u{idom} = obj.boundaryConditions.dirichletStep{idom}.reducedToFullVector(u{idom});
                intDOF = obj.interfaceDof(:,idom);
                R = R + w(idom)*obj.LHS{idom}(intDOF,:)*u{idom};
            end
        end
        
        function updateNeumanValues(obj,R)
            for idom=1:obj.nSubdomains(1)
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind = find(obj.boundaryConditions.neumannStep{idom}.neumann == obj.interfaceDof(idof,idom));
                    obj.boundaryConditions.neumannStep{idom}.neumann_values(ind) = R(idof);
                end
            end
        end

        function u = neumannStep(obj)
            for idom = 1:obj.nSubdomains(1)
                obj.computeForces(idom,obj.boundaryConditions.neumannStep)
                Kred    = obj.boundaryConditions.neumannStep{idom}.fullToReducedMatrix(obj.LHS{idom});
                Fred    = obj.boundaryConditions.neumannStep{idom}.fullToReducedVector(obj.RHS{idom});
                u{idom} = Kred\Fred;
            end
        end

        function uInt = updateDirichletStepBC(obj,u)
            uInt = obj.computeInterfaceDisp(u);
            obj.updateDirichletValues(uInt);
        end

        function uInt = computeInterfaceDisp(obj,u)
            w = [obj.theta,1- obj.theta];
             uInt=zeros(size(obj.interfaceDof,1),1);
             for idom=1:obj.nSubdomains(1)
                u{idom} = obj.boundaryConditions.neumannStep{idom}.reducedToFullVector(u{idom});
                intDOF = obj.interfaceDof(:,idom);
                uInt = uInt + w(idom)*u{idom}(intDOF);
            end
        end 

        function updateDirichletValues(obj,u)
            for idom=1:obj.nSubdomains(1)
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind = find(obj.boundaryConditions.dirichletStep{idom}.dirichlet == obj.interfaceDof(idof,idom));
                    obj.boundaryConditions.dirichletStep{idom}.dirichlet_values(ind) =...
                        obj.boundaryConditions.dirichletStep{idom}.dirichlet_values(ind) + 0.1*u(idof);
                end
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
