classdef AssemblerFun < handle

    properties (Access = private)
        fun
    end

    methods (Access = public)
        
        function obj = AssemblerFun(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, Aelem, test, trial)
            dofsTe = test.getConnec();
            dofsTr = trial.getConnec();            
            nElem = size(Aelem,3);
            nDofsTe = test.nDofs;
            nDofsTr = trial.nDofs;
            ndofElTe = test.nDofsElem;
            ndofElTr = trial.nDofsElem;

            res = zeros(ndofElTe*ndofElTr * nElem, 3);
            strt = 1;
            fnsh = nElem;
            for i = 1:ndofElTr
                dofsI = dofsTr(:,i);
                for j = 1:ndofElTe
                    dofsJ = dofsTe(:,j);
                    a = squeeze(Aelem(i,j,:));
                    matRes = [dofsI, dofsJ, a];
                    res(strt:fnsh,:) = matRes;
                    strt = strt + nElem;
                    fnsh = fnsh + nElem;
                end
            end
            A = sparse(res(:,1), res(:,2), res(:,3), nDofsTr, nDofsTe);

        end

        function V = assembleV(obj, F, fun)
            % Via indices
            dofConnec = obj.fun.getConnec();
            nDofsEl   = size(dofConnec,2);
            nDofs     = max(max(dofConnec)); %obj.fun.nDofs;
            nGaus     = size(F,2);
            nElem     = size(F,3);
            strt = 1;
            fnsh = nElem;
            res = zeros(nDofsEl * nElem, 2);
            for iDof = 1:nDofsEl
                for igaus = 1:nGaus
                    dofs = dofConnec(:,iDof);
                    c = squeeze(F(iDof,igaus,:));
                    matRes = [dofs, c];
                    res(strt:fnsh,:) = matRes;
                    strt = strt + nElem;
                    fnsh = fnsh + nElem;
                end
            end
            V = sparse(res(:,1), 1, res(:,2), nDofs, 1);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun    = cParams.fun;
        end

    end

end