classdef DomainDecompositionDofManager < handle

    properties (GetAccess = public, SetAccess = private)
        localGlobalDofConnec
        interfaceDof
        interfaceDom
        nDimf
    end

    properties (Access = private)
        nSubdomains
        meshDomain
        coarseMesh
        interfaceConnec
        locGlobConnec
        nBoundaryNodes
        nReferenceNodes
    end

    methods (Access = public)

        function obj = DomainDecompositionDofManager(cParams)
            obj.init(cParams)
            obj.createlocalGlobalDofConnec();
            obj.computeLocalInterfaceDof();
        end

        function [interfaceDof,interfaceDom] = getInterfaceDofDom(obj)
            interfaceDof = obj.interfaceDof;
            interfaceDom = obj.interfaceDom;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nSubdomains     = cParams.nSubdomains;
            obj.meshDomain      = cParams.meshDomain;
            obj.coarseMesh      = cParams.coarseMesh;
            obj.interfaceConnec = cParams.interfaceConnec;
            obj.locGlobConnec   = cParams.locGlobConnec;
            obj.nBoundaryNodes  = cParams.nBoundaryNodes;
            obj.nReferenceNodes = cParams.nReferenceNodes;
            obj.nDimf           = cParams.nDimf;
        end

        function  createlocalGlobalDofConnec(obj)
            ndimf = obj.nDimf;
            ndom  = obj.nSubdomains(1)*obj.nSubdomains(2);
            for dom = 1:ndom
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                localGlobalDofConnecDom = zeros(1,2);
                nodeG = obj.locGlobConnec(:,1,dom);
                nodeL = obj.locGlobConnec(:,2,dom);
                for iunkn = 1:ndimf
                    dofConec = [ndimf*(nodeG - 1) + iunkn ,  ndimf*(nodeL - 1) + iunkn] ;
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec ];
                end
                localGlobalDofConnec{row,col} = localGlobalDofConnecDom(2:end,:);
            end
            obj.localGlobalDofConnec = localGlobalDofConnec;
        end

        function computeLocalInterfaceDof(obj)
            intConec = reshape(obj.interfaceConnec',2,obj.nBoundaryNodes,[]);
            intConec = permute(intConec,[2 1 3]);
            nint = size(intConec,3);
            globaldof=0;
            ndimf = obj.nDimf;
            ndofs = obj.nReferenceNodes*ndimf;
            for iint=1:nint
                ndom = size(intConec,2); %length(intConec(1,:,iint));
                for idom = 1:ndom
                    dofaux=0;
                    nodesI = intConec(:,idom,iint);
                    dom = ceil(intConec(1,idom,iint)/obj.nReferenceNodes);
                    globaldof = (dom-1)*ndofs;
                    for iunkn=1:ndimf
                        DOF = ndimf*(nodesI - 1) + iunkn;
                        DOF = DOF-globaldof;
                        dofaux= [dofaux; DOF];
                    end
                    interfaceDof(:,idom,iint) = dofaux(2:end);
                    interfaceDom(iint,idom) = dom;
                    %                     globaldof = globaldof + (iint*(idom-1)+iint)*dim.ndofs;
                end
                %                 interfaceDof(:,iint) = dofaux(2:end);
                %                 globaldof = globaldof + dim.ndofs;
            end
            obj.interfaceDof = interfaceDof;
            obj.interfaceDom = interfaceDom;
        end

    end

end