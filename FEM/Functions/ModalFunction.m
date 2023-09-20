classdef ModalFunction < L2Function

    properties (Access = public)
        ndimf
    end

    properties (Access = private)
         fValues
%         mesh
        basisValues
        basisFunctions
        FEfun
        nbasis
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = ModalFunction(cParams)
            obj.init(cParams)
            obj.computeBasisFunctions();
        end

        function fxV = evaluate(obj, xGLoc)
            fxV = zeros(size(xGLoc));
            for ibasis = 1:obj.nbasis
                fI   = obj.fValues(ibasis);
                phiI = obj.basisFunctions{ibasis}.evaluate(xGLoc);
                fxV = fxV + phiI*fI;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.ndimf       = obj.mesh.ndim;
            obj.fValues     = cParams.fValues;
            obj.basisValues = cParams.basis;
            obj.nbasis  = numel(cParams.basis);
        end

%         function fvalue = dof2nodesFields(obj,basis)
%             nmodes = size(basis,2);
%             ndimf  = obj.ndimf;
%             nnode  = obj.mesh.nnodes;
%             %             for idimf=1:ndimf
%             % %                 aux=reshape(basis(idimf:ndimf:end,:),[nnode nmodes]);
%             %                aux=basis(idimf:ndimf:end,:);
%             %                fvalue(:,:,idimf)=aux;
%             %             end
% 
%             for imode=1:nmodes
%                 for idimf=1:ndimf
%                     fvalue(:,idimf,imode)=basis(idimf:ndimf:end,imode);
%                 end
%             end
%         end


        function computeBasisFunctions(obj)
            s.mesh  = obj.mesh;
            for ibasis=1:obj.nbasis
               s.fValues = obj.basisValues{ibasis};
               obj.basisFunctions{ibasis} = P1Function(s);
               obj.basisFunctions{ibasis}.plot
            end    
        end

    end

end