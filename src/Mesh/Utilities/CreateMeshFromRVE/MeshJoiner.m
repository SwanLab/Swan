classdef MeshJoiner < handle

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

        function obj = MeshJoiner(cParams)  
            obj.init(cParams)
        end

        function [mD,interfaceConnec,bdSB,localGlobalConnec,interfaceConnecReshaped,dMesh] = create(obj)
%             obj.createSubDomainMeshes();
            obj.createInterfaceSubDomainMeshes();
            obj.createDomainMesh();
            mD  = obj.meshDomain;
            interfaceConnec = obj.interfaceConnec;
            interfaceConnecReshaped = obj.interfaceConnecReshaped;
            bdSB = obj.interfaceMeshSubDomain;
            localGlobalConnec = obj.localGlobalConnec;
            dMesh = obj.domainMeshDisc;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nSubdomains = cParams.nsubdomains; %nx ny
            obj.meshReference = cParams.meshReference;
            obj.interfaceMeshReference  = obj.meshReference.createBoundaryMesh();
            obj.ninterfaces = length(obj.interfaceMeshReference);
            obj.tolSameNode = cParams.tolSameNode;
            obj.meshSubDomain = cParams.meshSbd;
        end

        function createInterfaceSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            for jDom = 1:nY
                for iDom = 1:nX
                    bIJ = obj.meshSubDomain{jDom,iDom}.createBoundaryMesh();
                    bD{jDom,iDom} = bIJ;
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
            coupling = InterfaceCoupling(s);
            coupling.compute();
            obj.interfaceConnec = coupling.interfaceConnec;
            
            
            s.interfaceConnec   = coupling.interfaceConnec;
            s.tolSameNode       = obj.tolSameNode;

            DMesh = DomainMeshComputer(s);
            DMesh.compute();
            obj.meshDomain        = DMesh.domainMesh;
            obj.localGlobalConnec = DMesh.localGlobalConnec;
            obj.domainMeshDisc    = DMesh.domainMeshDisc;

            obj.interfaceConnecReshaped = coupling.reshapeConecPerInterface();
            
        end
    end

end