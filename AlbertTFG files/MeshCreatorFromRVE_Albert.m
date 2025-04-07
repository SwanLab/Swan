classdef MeshCreatorFromRVE_Albert < handle

    properties (Access = public)

    end

    properties (Access = private)
        nSubdomains
        rSubdomains
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

        function obj = MeshCreatorFromRVE_Albert(cParams)  
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
            obj.nSubdomains   = cParams.nsubdomains; %nx ny
            obj.rSubdomains   = cParams.rsubdomains;
            obj.meshReference = cParams.meshReference;
            obj.tolSameNode   = cParams.tolSameNode;

            iMR = [];
            nin = [];

            for i = 1:size(obj.nSubdomains,2)
                imr = obj.meshReference(i).createBoundaryMesh();
                nin = cat(1, nin, length(imr) );
                iMR = cat(1, iMR, imr); 
            end

            obj.interfaceMeshReference = iMR;
            obj.ninterfaces            = nin;
            
        end

        function createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            %figure(2)
            for jDom = 1:nY
                for iDom = 1:nX
                    coordIJ = obj.computeSubdomainCoords(jDom,iDom);
                    mIJ     = obj.createSubdomainMesh(coordIJ,jDom,iDom);
                    %mIJ.plot();
                    mD{jDom,iDom} = mIJ;
                    % hold on
                end
            end
            obj.meshSubDomain = mD;
        end

        function concatenateMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            coord = [];
            connec = [];
            for jDom = 1:nX
                for iDom = 1:nY
                    m = obj.meshSubDomain{jDom,iDom};
                    coord  = [coord;m.coord];
                    connec = [connec;m.connec];
                end
            end
            s.coord  = coord;
            s.connec = connec;
            mAll = Mesh.create(s);

        end

        function L = computeReferenceMeshLength(obj,jDom,iDom)
            coord = obj.meshReference(jDom,iDom).coord;
            Lx = max(coord(:,1))-min(coord(:,1));
            Ly = max(coord(:,2))-min(coord(:,2));
            L = [Lx Ly];
        end

        function coord = computeSubdomainCoords(obj,jDom,iDom)
            coord0 = obj.meshReference(jDom, iDom).coord;
            L  = obj.computeReferenceMeshLength(jDom,iDom);
            Lx = L(1);
            Ly = L(2);
            coord(:,1) = coord0(:,1)+Lx*(iDom-1);
            coord(:,2) = coord0(:,2)+Ly*(jDom-1);
        end

        function m = createSubdomainMesh(obj,coord,jDom,iDom)
            connec0  = obj.meshReference(jDom,iDom).connec;
            s.coord  = coord;
            s.connec = connec0;
            m = Mesh.create(s);
        end

        function createInterfaceSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            
            for jDom = 1:nX
                for iDom = 1:nY
                    bIJ = obj.meshSubDomain{iDom,jDom}.createBoundaryMesh();
                    bD{jDom,iDom} = bIJ;
                  %   hold on
                  %   for iline=1:length(bIJ)
                  %       bIJ{iline}.mesh.plot();
                  %   end
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
            coupling = InterfaceCoupling_Albert(s);
            coupling.compute();
            obj.interfaceConnec = coupling.interfaceConnec;
            
            obj.interfaceConnecReshaped = coupling.interfaceConnecReshaped;
            
            s.interfaceConnec   = coupling.interfaceConnec;
            s.tolSameNode       = obj.tolSameNode;

            DMesh = DomainMeshComputer(s);
            DMesh.compute();
            obj.meshDomain        = DMesh.domainMesh;
            obj.localGlobalConnec = DMesh.localGlobalConnec;
            obj.domainMeshDisc    = DMesh.domainMeshDisc;
        end
    end

end