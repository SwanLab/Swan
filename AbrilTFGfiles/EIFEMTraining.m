classdef EIFEMTraining < handle

    properties (GetAccess = public, SetAccess = private)
        mesh
        Coarseorder
    end
    properties (Access = private)
        meshDomain
        DDdofManager
        domainIndices
        tolSameNode
        cellMeshes
        nSubdomains
        material
        geometryType
        levelSet
        unfittedMesh
        type
    end


    methods (Access = public)

        function obj = EIFEMTraining(cParams)
            obj.init(cParams);
        end


        function data=train(obj)
          
            obj.repeatMesh();  %create MeshDomain
            switch obj.type
                case 'continuous'
                    bMesh        = obj.meshDomain.createSingleBoundaryMesh();
                    s.lambdaFun  = LagrangianFunction.create(bMesh,obj.meshDomain.ndim,'P1');
                case 'discontinuous'
                    bMesh          = obj.meshDomain.createBoundaryMesh();
                    for i=1:numel(bMesh)
                        s.lambdaFun{i}=LagrangianFunction.create(bMesh{i}.mesh,obj.meshDomain.ndim, 'P1');
                    end
            end
            s.mesh          = obj.meshDomain;
            s.bMesh         = bMesh;
            s.uFun          = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            s.material      = obj.material;
            s.dirichletFun  = obj.createDirichletFunction(bMesh);
            s.type          = obj.type;
            e               = ElasticHarmonicExtension(s);
            [u,~,K,~]       = e.solve();

            [data.uSbd,data.LHSsbd] = obj.extractDomainData(u,K);

            % obj.print(u);
            % obj.print(data.uSbd);

             data.mesh= obj.mesh;
        end
    end

    
    methods (Access = private)

        function print(obj,T)
            for i=1:size(T,2)
                z.mesh      = obj.meshDomain;
                z.order     = 'P1';
                z.fValues   = reshape(T(:,i),z.mesh.ndim,[])';
                uFeFun = LagrangianFunction(z);%
                % uFeFun.print("Mesh3D"+num2str(i),'Paraview');
                uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(uFeFun);

                fvalues = [uMeshFun.innerMeshFunction.fValues;
                    uMeshFun.innerCutMeshFunction.fValues];

                s.coord = [uMeshFun.innerMeshFunction.mesh.coord;
                    uMeshFun.innerCutMeshFunction.mesh.coord];

                s.connec = [uMeshFun.innerMeshFunction.mesh.connec;
                    uMeshFun.innerCutMeshFunction.mesh.connec  + max(uMeshFun.innerMeshFunction.mesh.connec(:))];
                
                mh = Mesh.create(s);
                ss.mesh = mh;
                ss.fValues = fvalues;
                ss.order = 'P1';
                ss.ndimf = size(fvalues,2);
                u = LagrangianFunction(ss);
                u.print("LevelSet"+num2str(i),'Paraview')
            end

        end

        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.nSubdomains    = cParams.nSubdomains;
            obj.domainIndices  = cParams.domainIndices;
            obj.material       = cParams.material;
            % obj.levelSet       = cParams.levelSet;
            % obj.unfittedMesh   = cParams.unfittedMesh;
            obj.tolSameNode    = 1e-11;
            obj.Coarseorder    = cParams.Coarseorder;
            obj.type           = cParams.type;
        end

        function repeatMesh(obj)
            if sum(obj.nSubdomains > 1)>= 1
                bS  = obj.mesh.createBoundaryMesh();
                [mD,mSb,iC,lG,iCR,~] = obj.createMeshDomain(obj.mesh);
                obj.cellMeshes = mSb;
                obj.meshDomain = mD;
                obj.DDdofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,obj.mesh,iCR);
            else
                obj.cellMeshes= {obj.mesh};
                obj.meshDomain = obj.mesh;
            end
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            s.geometryType = obj.geometryType;
            m = MeshCreatorFromRVE.create(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function [u,lhs] = extractDomainData(obj,uC,LHS)
            if sum(obj.nSubdomains > 1)>= 1
                u   = obj.extractDomainDisplacements(uC);
                % lhs = obj.extractDomainLHS(LHS);
                lhs =1;
            else
                u = full(uC);
                lhs = LHS;
            end
        end

        function u = extractDomainDisplacements(obj,uC)
            ntest = size(uC,2);
            for i = 1:ntest
                if obj.mesh.ndim==2
                    uD    = obj.DDdofManager.global2local(uC(:,i));
                    ind   = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
                    u(:,i)= uD(:,ind);
                else
                    ind = obj.nSubdomains(1)*obj.nSubdomains(2)*(obj.domainIndices(3)-1)...
                        +obj.nSubdomains(1)*(obj.domainIndices(2)-1)+obj.domainIndices(1);
                    uD    = obj.DDdofManager.global2local(uC(:,i));
                    u(:,i)= uD(:,ind);
                end
            end
        end

        function lhs = extractDomainLHS(obj,LHS)
            lhs = obj.DDdofManager.global2localMatrix(LHS);
            if obj.mesh.ndim==2
                ind = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
            else
                ind = obj.nSubdomains(1)*obj.nSubdomains(2)*(obj.domainIndices(3)-1)...
                    +obj.nSubdomains(1)*(obj.domainIndices(2)-1)+obj.domainIndices(1);
            end
       
            lhs = lhs(:,:,ind);
        end

        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.meshDomain.nnodes;
            s.nDimf           = obj.meshDomain.ndim;

            if obj.mesh.ndim ==2
                d = DomainDecompositionDofManager(s);
            else
                d = DomainDecompositionDofManager3D(s);
            end
        end

        function dF = createDirichletFunction(obj,bMesh)
            s.mesh = bMesh;
            s.order= obj.Coarseorder;
            s.type = obj.type;
            cf = CoarseFunctions(s);
            dF = cf.getAnalytical();
        end

    end


    
end
