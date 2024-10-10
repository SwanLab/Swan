classdef DomainDecompositionDofManager < handle

    properties (Access = private)
        nDof
        nReferenceDof
        localGlobalDof
        interfaceDof
        interfaceDom
    end

    properties (Access = private)
        nSubdomains
        interfaceConnec
        locGlobConnec
        nBoundaryNodes
        nReferenceNodes
        nDimf
        nNodes
    end

    methods (Access = public)

        function init(obj,cParams)
            obj.nSubdomains     = cParams.nSubdomains;
            obj.interfaceConnec = cParams.interfaceConnec;
            obj.locGlobConnec   = cParams.locGlobConnec;
            obj.nBoundaryNodes  = cParams.nBoundaryNodes;
            obj.nReferenceNodes = cParams.nReferenceNodes;
            obj.nNodes          = cParams.nNodes;
            obj.nDimf           = cParams.nDimf;
            obj.nDof            = obj.nNodes*obj.nDimf;
            obj.nReferenceDof   = obj.nReferenceNodes*obj.nDimf;
        end        

        function obj = DomainDecompositionDofManager(cParams)
            obj.init(cParams)
            obj.createlocalGlobalDofConnec();
            obj.computeLocalInterfaceDof();
        end

        function f = scaleInterfaceValues(obj,f,w)
            nint = size(obj.interfaceDof,3);
            weight = [w,1-w];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    dof = obj.interfaceDof(:,idom,iint);
                    f(dof,dom) = weight(idom)* f(dof,dom);
                end
            end
        end

        function fG = local2global(obj,fL)
            fG   = zeros(obj.nDof,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    lDof  = obj.localGlobalDof{jdom,idom}(:,1);
                    gDof = obj.localGlobalDof{jdom,idom}(:,2);
                    fG(lDof,ind) = fL(gDof,ind);
                    ind=ind+1;
                end
            end
        end

        function fL = global2local(obj,fG)
            ndimf  = obj.nDimf;
            fL     = zeros(obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    lGconnec = obj.localGlobalDof{jdom,idom};
                    fL(lGconnec(:,2),ind) = fG(lGconnec(:,1));
                    ind=ind+1;
                end
            end
        end
    end

    methods (Access = private)

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
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec];
                end
                localGlobalDofConnec{row,col} = localGlobalDofConnecDom(2:end,:);
            end
            obj.localGlobalDof = localGlobalDofConnec;
        end

        function computeLocalInterfaceDof(obj)
            intConec = reshape(obj.interfaceConnec',2,obj.nBoundaryNodes,[]);
            intConec = permute(intConec,[2 1 3]);
            nint = size(intConec,3);
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
                end
            end
            obj.interfaceDof = interfaceDof;
            obj.interfaceDom = interfaceDom;
        end

    end

end