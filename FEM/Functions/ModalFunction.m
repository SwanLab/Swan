classdef ModalFunction < L2Function

    properties (Access = public)
        ndimf
        nbasis
        fValues
        basisFunctions
    end

    properties (Access = private)
        basisValues
        functionType
    end


    methods (Access = public)

        function obj = ModalFunction(cParams)
            obj.init(cParams)
            obj.computeBasisFunctions();
        end

        function fxV = evaluate(obj, xGLoc)
            nelem=obj.mesh.nelem;
            sizeaux= [obj.ndimf,size(xGLoc,2),nelem];
            fxV = zeros(sizeaux);
            for ibasis = 1:obj.nbasis
                fI   = obj.fValues(ibasis);
                phiI = obj.basisFunctions{ibasis}.evaluate(xGLoc);
                fxV = fxV + phiI*fI;
            end
        end

        function fxV = evaluateBasisFunctions(obj,xGLoc)
            for ibasis=1:obj.nbasis
               phiI = obj.basisFunctions{ibasis}.evaluate(xGLoc);
               fxV{ibasis}=phiI;
            end    
        end

         function plot(obj)
            p1DiscFun = obj.project('P1D');
            p1DiscFun.plot();
        end

    end

    methods (Access = public, Static)

        function MF = create(mesh,basis,functionType)
            nbasis         = numel(basis);
            s.fValues      = zeros(nbasis,1);
            s.mesh         = mesh;
            s.basis        = basis;
            s.functionType = functionType;
            MF = ModalFunction(s);
        end

        function fS = times(f1,f2)
            fS = f1.fValues.*f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = ModalFunction(s);
        end

        function fS = sum(f1,f2)
            fS = f1.fValues+f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = ModalFunction(s);
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.ndimf        = obj.mesh.ndim;
            obj.fValues      = cParams.fValues;
            obj.basisValues  = cParams.basis;
            obj.functionType = cParams.functionType;
            obj.nbasis       = numel(cParams.basis);
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
               s.fValues      = obj.basisValues{ibasis};
               s.functionType = obj.functionType{ibasis};
               obj.basisFunctions{ibasis} = FunctionFactory.create(s);
%                obj.basisFunctions{ibasis} = P1Function(s);
%                obj.basisFunctions{ibasis}.plot
            end    
        end

    end

end