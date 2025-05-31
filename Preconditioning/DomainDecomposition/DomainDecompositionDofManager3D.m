classdef DomainDecompositionDofManager3D < handle

    properties (Access = public)
        interfaceConnec
        interfaceConnecReshaped
        interfaceDom
        intConecLocal
    end

    properties (Access = private)
        nDof
        nReferenceDof
        localGlobalDof
        interfaceDof        
    end

    properties (Access = private)
        nSubdomains
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
            obj.interfaceConnecReshaped = cParams.interfaceConnecReshaped;
            obj.locGlobConnec   = cParams.locGlobConnec;
            obj.nBoundaryNodes  = cParams.nBoundaryNodes;
            obj.nReferenceNodes = cParams.nReferenceNodes;
            obj.nNodes          = cParams.nNodes;
            obj.nDimf           = cParams.nDimf;
            obj.nDof            = obj.nNodes*obj.nDimf;
            obj.nReferenceDof   = obj.nReferenceNodes*obj.nDimf;
        end        

        function obj = DomainDecompositionDofManager3D(cParams)
            obj.init(cParams)
            obj.createlocalGlobalDofConnec();
            obj.computeLocalInterfaceDof();
        end

        function f = scaleInterfaceValues(obj,f,w)
            nint = numel(obj.interfaceDof);
            weight = [w,1-w];
            for iint = 1:nint
                dofI = obj.interfaceDof{iint};
                ndom = size(dofI,2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    dof = dofI(:,idom);
                    f(dof,dom) = weight(idom)* f(dof,dom);
                end
            end
        end

        function m = scaleInterfaceValuesMatrix(obj,m,w)
            nint = numel(obj.interfaceDof);
            weight = [w,1-w];
            for iint = 1:nint
                dofI = obj.interfaceDof{iint};
                ndom = size(dofI,2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    dof = dofI(:,idom);
                    m(dof,dof,dom) = weight(idom)* m(dof,dof,dom);
                end
            end
        end

        function fG = local2global(obj,fL)
            fG   = zeros(obj.nDof,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for kdom = 1:obj.nSubdomains(3)
                for jdom = 1: obj.nSubdomains(2)
                    for idom = 1: obj.nSubdomains(1)
                        lDof  = obj.localGlobalDof{jdom,idom,kdom}(:,1);
                        gDof = obj.localGlobalDof{jdom,idom,kdom}(:,2);
                        fG(lDof,ind) = fL(gDof,ind);
                        ind=ind+1;
                    end
                end
            end
        end

        function fL = global2local(obj,fG)
            ndimf  = obj.nDimf;
            fL     = zeros(obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for kdom=1:obj.nSubdomains(3)
                for jdom = 1:obj.nSubdomains(2)
                    for idom = 1:obj.nSubdomains(1)
                        lGconnec = obj.localGlobalDof{jdom,idom,kdom};
                        fL(lGconnec(:,2),ind) = fG(lGconnec(:,1));
                        ind=ind+1;
                    end
                end
            end
        end

        function mL = global2localMatrix(obj,mG)
            ndimf  = obj.nDimf;
            mL     = zeros(obj.nReferenceNodes*ndimf,obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for kdom = 1:obj.nSubdomains(3)
                for jdom = 1:obj.nSubdomains(2)
                    for idom = 1:obj.nSubdomains(1)
                        lGconnec = obj.localGlobalDof{jdom,idom,kdom};
                        mL(lGconnec(:,2),lGconnec(:,2),ind) = mG(lGconnec(:,1),lGconnec(:,1));
                        ind=ind+1;
                    end
                end
            end
        end
        
    end

    methods (Access = private)

     function  createlocalGlobalDofConnec(obj)
            ndimf = obj.nDimf;
            ndom  = obj.nSubdomains(1)*obj.nSubdomains(2)*obj.nSubdomains(3);
            for dom = 1:ndom
                kInd   = ceil(dom/(obj.nSubdomains(1)*obj.nSubdomains(2)));
                pageD  = dom-(kInd-1)*obj.nSubdomains(1)*obj.nSubdomains(2);
                row = ceil(pageD/obj.nSubdomains(1));
                col = pageD-(row-1)*obj.nSubdomains(1);
%                 row = ceil(dom/obj.nSubdomains(1));
%                 col = dom-(row-1)*obj.nSubdomains(1);
                localGlobalDofConnecDom = zeros(1,2);
                nodeG = obj.locGlobConnec(:,1,dom);
                nodeL = obj.locGlobConnec(:,2,dom);
                for iunkn = 1:ndimf
                    dofConec = [ndimf*(nodeG - 1) + iunkn ,  ndimf*(nodeL - 1) + iunkn] ;
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec];
                end
                localGlobalDofConnec{row,col,kInd} = localGlobalDofConnecDom(2:end,:);
            end
            obj.localGlobalDof = localGlobalDofConnec;
        end

        function computeLocalInterfaceDof(obj)
                       
            intConecResh = obj.interfaceConnecReshaped;



            nint = numel(intConecResh);
            ndimf = obj.nDimf;
            ndofs = obj.nReferenceNodes*ndimf;
            for iint=1:nint
                intConecIint = intConecResh{iint};                
                ndom = size(intConecIint,2); %length(intConec(1,:,iint));
                intConecL = zeros(size(intConecIint));
                interfaceDof = [];
                for idom = 1:ndom
                    dofaux=0;
                    
                    nodesI = intConecIint(:,idom);
                    dom = ceil(intConecIint(1,idom)/obj.nReferenceNodes);
                    globaldof = (dom-1)*ndofs;
                    for iunkn=1:ndimf
                        DOF = ndimf*(nodesI - 1) + iunkn;
                        DOF = DOF-globaldof;
                        dofaux= [dofaux; DOF];
                    end
                    interfaceDof(:,idom) = dofaux(2:end);
                    interfaceDom(iint,idom) = dom;
                    intConecL(:,idom) = nodesI - (dom-1)*obj.nReferenceNodes;
                end
                intConecLiint{iint} = intConecL(:,idom);
                interfaceDofIint{iint} = interfaceDof;
            end
            obj.interfaceDof  = interfaceDofIint;
            obj.interfaceDom  = interfaceDom;
            obj.intConecLocal = intConecLiint;

        end

    end

end