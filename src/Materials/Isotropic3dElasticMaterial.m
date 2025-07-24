classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = public)

        function C = evaluate(obj,xV)
            [l,m] = computeLameParameters(obj);
            lambda = l.evaluate(xV);
            mu = m.evaluate(xV);

            N = obj.ndim;
            nGauss = size(mu,2);
            nElem  = m.mesh.nelem;
            lambda = reshape(lambda,[1 1 1 1 nGauss nElem]);
            mu     = reshape(mu,[1 1 1 1 nGauss nElem]);
            I      = repmat(eye4D(N),[1 1 1 1 nGauss nElem]);
            IxI    = repmat(kronEye(N),[1 1 1 1 nGauss nElem]);
            C = 2*mu.*I + lambda.*IxI;
        end
        
    end
    
end

