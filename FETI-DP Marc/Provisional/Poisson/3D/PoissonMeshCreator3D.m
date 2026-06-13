classdef PoissonMeshCreator3D < handle

    properties (Access = private)
        nSubdomains
        meshReference
        interfaceMeshReference
        ninterfaces
        tolSameNode
    end

    properties (Access = private)
        meshSubDomain
        interfaceMeshSubDomain
        meshDomain
        interfaceConnec
        localGlobalConnec
        interfaceConnecReshaped
        domainMeshDisc
    end

    methods (Access = public)

        function obj = PoissonMeshCreator3D(cParams)  
            obj.init(cParams)
        end

        function [mD,mSD,interfaceConnec,bdSB,localGlobalConnec,interfaceConnecReshaped,dMesh] = create(obj)
            obj.createSubDomainMeshes();
            obj.createInterfaceSubDomainMeshes();
            obj.createDomainMesh();
            mD  = obj.meshDomain;
            mSD = obj.meshSubDomain;
            interfaceConnec = obj.interfaceConnec;
            interfaceConnecReshaped = obj.interfaceConnecReshaped;
            bdSB = obj.interfaceMeshSubDomain;
            localGlobalConnec = obj.localGlobalConnec;
            dMesh = obj.domainMeshDisc;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nSubdomains = cParams.nsubdomains;
            obj.meshReference = cParams.meshReference;
            obj.interfaceMeshReference  = obj.meshReference.createBoundaryMesh();
            obj.ninterfaces = length(obj.interfaceMeshReference);
            obj.tolSameNode = cParams.tolSameNode;
        end

        function createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nZ = obj.nSubdomains(3);
            for kDom = 1: nZ
                for jDom = 1:nY
                    for iDom = 1:nX
                        coordIJ = obj.computeSubdomainCoords(jDom,iDom,kDom);
                        mIJK    = obj.createSubdomainMesh(coordIJ);
                        mD{jDom,iDom,kDom} = mIJK;
                    end
                end
            end
            obj.meshSubDomain = mD;
        end

        function L = computeReferenceMeshLength(obj)
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1))-min(coord(:,1));
            Ly = max(coord(:,2))-min(coord(:,2));
            Lz = max(coord(:,3))-min(coord(:,3));
            L = [Lx Ly Lz];
        end

        function coord = computeSubdomainCoords(obj,jDom,iDom,kDom)
            coord0 = obj.meshReference.coord;
            L  = obj.computeReferenceMeshLength();
            Lx = L(1);
            Ly = L(2);
            Lz = L(3);
            coord(:,1) = coord0(:,1)+Lx*(iDom-1);
            coord(:,2) = coord0(:,2)+Ly*(jDom-1);
            coord(:,3) = coord0(:,3)+Lz*(kDom-1);
        end

        function m = createSubdomainMesh(obj,coord)
            connec0  = obj.meshReference.connec;
            s.coord  = coord;
            s.connec = connec0;
            m = Mesh.create(s);
        end

        function createInterfaceSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nZ = obj.nSubdomains(3);
            for kDom = 1:nZ
                for jDom = 1:nY
                    for iDom = 1:nX
                        bIJ = obj.meshSubDomain{jDom,iDom,kDom}.createBoundaryMesh();
                        bD{jDom,iDom,kDom} = bIJ;
                    end
                end
            end
            obj.interfaceMeshSubDomain = bD;
        end

        function createDomainMesh(obj)
            s.nSubdomains   = obj.nSubdomains;
            s.meshReference = obj.meshReference;
            s.interfaceMeshSubDomain = obj.interfaceMeshSubDomain;
            s.ninterfaces   = obj.ninterfaces;
            s.meshSubDomain = obj.meshSubDomain;
            s.tolSameNode = obj.tolSameNode;
            coupling = InterfaceCoupling3D(s);
            coupling.compute();
            obj.interfaceConnec = coupling.interfaceConnec;
            s.interfaceConnec   = coupling.interfaceConnec;
            s.tolSameNode       = obj.tolSameNode;
            DMesh = DomainMeshComputer3D(s);
            DMesh.compute();
            obj.meshDomain        = DMesh.domainMesh;
            obj.localGlobalConnec = DMesh.localGlobalConnec;
            obj.domainMeshDisc    = DMesh.domainMeshDisc;
            try
                obj.interfaceConnecReshaped = coupling.reshapeConecPerInterface();
            catch
                obj.interfaceConnecReshaped = {};
            end
        end
    end

end
