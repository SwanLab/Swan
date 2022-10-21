classdef Gradient < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fun
        mesh
        quadrature
    end
    
    methods (Access = public)
        
        function obj = Gradient()
            obj.init()
        end

        function precompute(obj)
            % fun.fValues still unknown
        end

        function gradFun = compute(obj, fun, quad, mesh)
            % fun.fValues known
            dNdx  = fun.computeCartesianDerivatives(quad, mesh);
            nDimf = fun.ndimf;
            nDims = size(dNdx, 1); % derivX, derivY (mesh-related?)
            nNode = size(dNdx, 2);
            nElem = size(dNdx, 3);
            nGaus = size(dNdx, 4);
            
            grad = zeros(nDims,nDimf, nElem, nGaus);
            for iGaus = 1:nGaus
                dNdx_g = dNdx(:,:,:,iGaus);
                for iDims = 1:nDims
                    for iNode = 1:nNode
                        dNdx_i = squeeze(dNdx_g(iDims, iNode,:));
                        nodes = mesh.connec(:,iNode);
                        f = fun.fValues(nodes,:);
                        p = (dNdx_i.*f)';
                        pp(1,:,:) = p;
                        grad(iDims,:,:,iGaus) = grad(iDims,:,:,iGaus) + pp;
                    end
                end
            end
            fVR = reshape(grad, [nDims*nDimf,nElem, nGaus]);
            s.fValues = permute(fVR, [1 3 2]);
%             s.ndimf      = nDimf;
            s.quadrature = quad;
            gradFun = FGaussDiscontinuousFunction(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
%             obj.fun        = cParams.fun;
%             obj.mesh       = cParams.mesh;
%             obj.quadrature = cParams.quadrature;
        end
        
    end
    
end