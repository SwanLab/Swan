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
            dofs   = obj.fun.computeDofConnectivity();
            nDofsE = size(dofs,1);
            nDofs  = numel(obj.fun.fValues);
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
        end

    end

end