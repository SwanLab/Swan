classdef AssemblerFun < handle

    properties (Access = private)
        fun
        connec
    end

    methods (Access = public)
        
        function obj = AssemblerFun(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, Aelem)
            dofs = obj.computeFunctionDofs();
            nDimf  = obj.fun.ndimf;
            nNodes = size(obj.fun.fValues,1);
            nNodeE = size(obj.connec,2);
            nDofsE = nNodeE*nDimf;
            nDofs  = nNodes*nDimf;
            A = sparse(nDofs,nDofs);
            for i = 1:nDofsE
                for j = 1:nDofsE
                    a = squeeze(Aelem(i,j,:));
                    A = A + sparse(dofs(i,:),dofs(j,:),a,nDofs,nDofs);
                end
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun    = cParams.fun;
            obj.connec = cParams.connec;
        end
        
        function dofConnec = computeFunctionDofs(obj)
            conne  = obj.connec;
            nDimf  = obj.fun.ndimf;
            nNode  = size(conne, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(conne,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = conne(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

    end

end