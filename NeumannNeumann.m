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
        Fext
        weight
        interfaceNeumanDof
        interfaceDirichletDof
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
            obj.interfaceDof = obj.computeLocalInterfaceDof();
            obj.computeSubdomainLHS();
            
            obj.createSubdomainBoundaryConditions();

            obj.interfaceNeumanDof = obj.identifyNeumanInterfaceDof();
            obj.interfaceDirichletDof = obj.identifyDirichletInterfaceDof();

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
            obj.weight   = 0.5;
            obj.theta    = 0.1;
        end

        function createReferenceMesh(obj)
            %             filename   = 'lattice_ex1';
            %             a.fileName = filename;
            %             femD       = FemDataContainer(a);
            %             mS         = femD.mesh;
            %             bS         = mS.createBoundaryMesh();
            % Generate coordinates
            x1 = linspace(0,1,12);
            x2 = linspace(0,0.5,12);
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

%             newmanBc{2}=[];
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
            dim = obj.getFunDims(obj.displacementFun{i});
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

        function dim = getFunDims(obj,disp)
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
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

        function [Fext,RHS] = computeForces(obj,boundaryConditions,material,mesh,disp,LHS)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims(disp);
            s.BC       = boundaryConditions;
            s.mesh     = mesh;
            s.material = material;
            %             s.globalConnec = obj.displacementFun{i}.connec;
            s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            Fext = rhs;
            R = RHSint.computeReactions(LHS);
            %             obj.variables.fext = rhs + R;
            RHS = rhs+R;
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
                dim = obj.getFunDims(obj.displacementFun{i});
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

        function ind = identifyNeumanInterfaceDof(obj)
            for idom=1:obj.nSubdomains(1)
                bN = obj.boundaryConditions.neumannStep{idom};
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind(idof,idom) = find(bN.neumann == obj.interfaceDof(idof,idom));
                end
            end
        end

        function ind = identifyDirichletInterfaceDof(obj)
             for idom=1:obj.nSubdomains(1)
                bD = obj.boundaryConditions.dirichletStep{idom};
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind(idof,idom) = find(bD.dirichlet == obj.interfaceDof(idof,idom));
                end
            end
        end

     
        function NeumannNeumannSolver(obj)
            tol=1e-8;
            e=1;
            iter=1;
            uInt = zeros(size(obj.interfaceDof,1),1); 
            while e(iter)>tol

                [uD,RHS] = obj.solveFEM(obj.boundaryConditions.dirichletStep);
                R = obj.computeInterfaceResidual(uD,RHS);
                obj.updateNeumanValues(R);
                
                [uN,~] = obj.solveFEM(obj.boundaryConditions.neumannStep);
                uInt_new = obj.updateDirichletStepBC(uN);

                iter=iter+1;
                e(iter) = norm(uInt_new-uInt);
                uInt = uInt_new;

%                 [uD,RHS] = obj.dirichletStep();
                
               % for idom = 1:obj.nSubdomains(1)
              %      uplot = reshape(obj.displacementFun{idom}.fValues',1,[])';
%                     uplot = obj.boundaryConditions.neumannStep{idom}.reducedToFullVector(uN{idom});
              %      obj.plotSolution(uplot,obj.meshSubDomain{idom}, idom,iter)
              %  end
%                 iter=iter+1;
                
%                 uN = obj.neumannStep();
                
                
              %  for idom = 1:obj.nSubdomains(1)
%                     uplot = obj.boundaryConditions.neumannStep{idom}.reducedToFullVector(uN{idom});
              %      uplot = uN{idom};
              %      obj.plotSolution(uplot,obj.meshSubDomain{idom}, idom,iter)
              %  end
               
%                 iter=iter+1;
%                 e(iter) = norm(uInt_new-uInt);
%                 uInt = uInt_new;

%                iter=iter+1;
                %                 end
            end
        end

        function [uD,RHS] = dirichletStep(obj)
            for idom = 1:obj.nSubdomains(1)
                bD = obj.boundaryConditions.dirichletStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                disp = obj.displacementFun{idom};
                mat  = obj.material{idom};
                [Fext,RHS{idom}] = obj.computeForces(bD,mat,mesh,disp,lhs);
                Kred    = bD.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bD.fullToReducedVector(RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bD.reducedToFullVector(u{idom});
                uD{idom} = reshape(u{idom},2,[])';                
            end
        end
        
        function u = neumannStep(obj)
            for idom = 1:obj.nSubdomains(1)
                bN = obj.boundaryConditions.neumannStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                mat  = obj.material{idom};                
                disp = obj.displacementFun{idom};
                [obj.Fext,obj.RHS{idom}] = obj.computeForces(bN,mat,mesh,disp,lhs);
                Kred    = bN.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bN.fullToReducedVector(obj.RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bN.reducedToFullVector(u{idom});
            end
        end

        function [u,RHS] = solveFEM(obj,bc)
            for idom = 1:obj.nSubdomains(1)
%                 bc = obj.boundaryConditions.neumannStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                mat  = obj.material{idom};                
                disp = obj.displacementFun{idom};
                bc_dom = bc{idom};
                [obj.Fext,RHS{idom}] = obj.computeForces(bc_dom,mat,mesh,disp,lhs);
                Kred    = bc_dom.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bc_dom.fullToReducedVector(RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bc_dom.reducedToFullVector(u{idom});
                u{idom} = reshape(u{idom},2,[])'; 
            end
        end
        
        function R2 = computeInterfaceResidual(obj,u,RHS)   
            
            w = [obj.weight,1-obj.weight];
%              w  = [1,1];
             R=zeros(size(obj.interfaceDof,1),1);
             R2=zeros(size(obj.interfaceDof,1),1);
             for idom=1:obj.nSubdomains(1)
                unodal = reshape(u{idom}',1,[])'; 
                interfaceDOF = obj.interfaceDof(:,idom);
                 R2 = R2 + w(idom)*(RHS{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal ...
                                  + obj.LHS{idom}(interfaceDOF,interfaceDOF)*unodal(interfaceDOF));
%                R = R + w(idom)*(obj.Fext{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal);
%                 aaa = norm(R2-R);
            end
        end
        
        function updateNeumanValues(obj,R)
            for idom=1:obj.nSubdomains(1)
                ind = obj.interfaceNeumanDof(:,idom);
                obj.boundaryConditions.neumannStep{idom}.neumann_values(ind) = R;
            end
        end


        function uIntNew = updateDirichletStepBC(obj,u)
            uInt = obj.computeInterfaceDisp(u);
            uIntNew = obj.updateDirichletValues(uInt);
        end

        function uInt = computeInterfaceDisp(obj,u)
            w = [obj.weight,1- obj.weight];
%             w  = [1,1];
             uInt=zeros(size(obj.interfaceDof,1),1);
             for idom=1:obj.nSubdomains(1)
                unodal = reshape(u{idom}',1,[])';
                interfaceDOF = obj.interfaceDof(:,idom);
                uInt = uInt + w(idom)*unodal(interfaceDOF);
            end
        end 

        function uIntNew = updateDirichletValues(obj,u)
            for idom=1:obj.nSubdomains(1)
                bD{idom} = obj.boundaryConditions.dirichletStep{idom};
                ind = obj.interfaceDirichletDof(:,idom);
                bD{idom}.dirichlet_values(ind) = bD{idom}.dirichlet_values(ind) + obj.theta*u;
           
            end
            obj.boundaryConditions.dirichletStep = bD;
            uIntNew = bD{end}.dirichlet_values(obj.interfaceDirichletDof(:,end));
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
