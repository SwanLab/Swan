classdef AssemblerFun < handle

    properties (Access = private)
        fun
    end

    methods (Access = public)
        
        function obj = AssemblerFun(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, Aelem, f1, f2)
            dofsF1 = f1.getDofConnec();
            if isequal(f1, f2)
                dofsF2 = dofsF1;
            else
                dofsF2 = f2.getDofConnec();
            end
            nDofs1     = numel(f1.fValues);
            nDofs2     = numel(f2.fValues);
            ndofsElem1 = size(Aelem, 1);
            ndofsElem2 = size(Aelem, 2);
            
            [iElem, jElem] = meshgrid(1:ndofsElem1, 1:ndofsElem2);
            iElem = iElem(:);
            jElem = jElem(:);

            dofsI = dofsF1(:, iElem);
            dofsJ = dofsF2(:, jElem);

            rowIdx = dofsI(:);
            colIdx = dofsJ(:);
            Aval   = permute(Aelem,[3 1 2]);
            values = Aval(:);

            A = sparse(rowIdx, colIdx, values, nDofs1, nDofs2);
        end

        function V = assembleV(obj, F, fun)
            % Via indices
            dofConnec = obj.fun.getDofConnec();
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