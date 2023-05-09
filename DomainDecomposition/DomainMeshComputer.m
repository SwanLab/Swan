classdef DomainMeshComputer < handle

    properties (Access = public)
        domainMesh
    end

    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
        meshSubDomain
        interfaceConnec
    end

    properties (Access = private)
        connecGlob
        coordGlob
    end

    methods (Access = public)

        function obj = DomainMeshComputer(cParams)
            obj.init(cParams)

        end

        function compute(obj)
            obj.createGlobalConnec();
            obj.createGlobalCoord();
            obj.updateGlobalConnec();
            obj.updateGlobalCoord();
            s.coord=obj.coordGlob;
            s.connec=obj.connecGlob;
            obj.domainMesh=Mesh(s); 
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.meshReference   = cParams.meshReference;
            obj.nSubdomains     = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces     = cParams.ninterfaces;
            obj.meshSubDomain   = cParams.meshSubDomain;
            obj.interfaceConnec = cParams.interfaceConnec;
        end

        function createGlobalConnec(obj)
            % we simply create a connectivity matrix as if no nodes are
            % shared, just summing nnodes. We need that to create the
            % global connectivity matrix that considers shared nodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nelem = obj.meshReference.nelem;
            nnodeElem=obj.meshReference.nnodeElem;
            nnodes  = obj.meshReference.nnodes;
            connec0 =obj.meshReference.connec;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            connecGlob=zeros(nX*nY*nelem,nnodeElem);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear= nX*(jDom-1)+iDom;
                    rowIn=(indLinear-1)*nelem+1;
                    rowEnd=indLinear*nelem;
                    connecGlob(rowIn:rowEnd,:)=connec0+nnodes*(indLinear-1);
                end
            end
            obj.connecGlob=connecGlob;
        end

        function  createGlobalCoord(obj)
            % we simply create the coordinate matrix as if no nodes are
            % shared, just summing nnodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            ndim     = obj.meshReference.ndim;
            nnodes  = obj.meshReference.nnodes;
            meshsd = obj.meshSubDomain;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            coordGlob=zeros(nX*nY*nnodes,ndim);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear= nX*(jDom-1)+iDom;
                    rowIn=(indLinear-1)*nnodes+1;
                    rowEnd=indLinear*nnodes;
                    coordGlob(rowIn:rowEnd,:)=meshsd{jDom,iDom}.coord;
                end
            end
            obj.coordGlob=coordGlob;
        end

        function   updateGlobalConnec(obj)
            updtConnecGlob = obj.connecGlob;
            rCconnec       = obj.interfaceConnec;
            nref           = size(rCconnec,1);
            ncopies        = size(rCconnec,2);
            for iref=1:nref
                for icopy=2:ncopies
                    updtConnecGlob(updtConnecGlob==rCconnec(iref,icopy))=rCconnec(iref,1);
                    updtConnecGlob(updtConnecGlob>rCconnec(iref,icopy))=updtConnecGlob(updtConnecGlob>rCconnec(iref,icopy))-1;
                    rCconnec(rCconnec>rCconnec(iref,icopy))=rCconnec(rCconnec>rCconnec(iref,icopy))-1;
                end
            end
            obj.connecGlob=updtConnecGlob;
        end

        function  updateGlobalCoord(obj)
            obj.coordGlob  = unique(obj.coordGlob,'rows','stable');
        end
    end
end