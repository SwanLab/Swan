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
    end


    methods (Access = public)

        function obj = EIFEMTraining(cParams)
            obj.init(cParams);
        end


        function data=train(obj)
            obj.repeatMesh();  %create MeshDomain
            bMesh = obj.meshDomain.createSingleBoundaryMesh();
            s.mesh=bMesh;
            s.type='continuous';
            cf=CoarseFunctions(s);
            f=cf.getAnalytical();
            dirchletFun = f;

            s.mesh         = obj.meshDomain;
            s.uFun         = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            s.lambdaFun    = LagrangianFunction.create(bMesh,obj.meshDomain.ndim,'P1');
            s.material     = obj.material;
            s.dirichletFun = dirchletFun;
            e  = ElasticHarmonicExtension(s);
            [u,~,K,~] = e.solve();

            [data.uSbd,data.LHSsbd] = obj.extractDomainData(u,K);

            obj.print(data.uSbd);

             data.mesh= obj.mesh;
             data.Coarseorder= obj.Coarseorder;
        end
    end

    
    methods (Access = private)

        function print(obj,T)
           z.mesh      = obj.mesh;
           z.order     = 'P1';
           z.fValues   = reshape(T(:,1),[obj.mesh.ndim,obj.mesh.nnodes])';
            uFeFun = LagrangianFunction(z);%            
            
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
            ss.ndimf = size(fvalues,2)

            u = LagrangianFunction(ss);

            u.print(["dsafaf"],'Paraview')

        end

        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.nSubdomains    = cParams.nSubdomains;
            obj.domainIndices  = cParams.domainIndices;
            obj.material       = cParams.material;
            obj.levelSet       = cParams.levelSet;
            obj.tolSameNode    = 1e-10;
            obj.Coarseorder    = 1;
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
                lhs = obj.extractDomainLHS(LHS);
            else
                u = full(uC);
                lhs = LHS;
            end
        end

        function u = extractDomainDisplacements(obj,uC)
            ntest = size(uC,2);
            for i = 1:ntest
                uD    = obj.DDdofManager.global2local(uC(:,i));
                ind   = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
                u(:,i)= uD(:,ind);
            end
        end

        function lhs = extractDomainLHS(obj,LHS)
            lhs = obj.DDdofManager.global2localMatrix(LHS);
            ind = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
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
            d = DomainDecompositionDofManager(s);
        end

    end


end
