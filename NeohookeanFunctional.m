classdef NeohookeanFunctional < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        lambda
        mu
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = NeohookeanFunctional(cParams)
            obj.init(cParams)
            
        end

        function val = compute(obj, uFun)
            quad = Quadrature.create(obj.mesh, 2);
            xG = quad.posgp;

            nPoints  = quad.ngaus;
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;
            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);

            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1/nPoints;
            I33(2,2,:,:) = 1/nPoints;
            I33(3,3,:,:) = 1/nPoints;

            F = I33 + GradU; % deformation gradient
            Ft = permute(F, [2 1 3 4]);
            
            C = pagemtimes(Ft,F);
            trC = C(1,1,:,:) +C(2,2,:,:)+C(3,3,:,:);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);

            val = obj.mu/2*(trC - 3) - obj.mu*log(jac) + obj.lambda/2*(jac-1).^2;
            val = squeeze(sum(val,3));
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.mesh   = cParams.mesh;
        end
        
    end
    
end