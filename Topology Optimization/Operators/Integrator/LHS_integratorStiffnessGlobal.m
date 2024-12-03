classdef LHS_integratorStiffnessGlobal < handle

    properties (Access = public)

    end

    properties (Access = private)
        material
        test
        trial
        quadrature
        quadratureOrder
        mesh
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = LHS_integratorStiffnessGlobal(cParams)
            obj.init(cParams)
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            quad  = obj.quadrature;
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            Cmat  = obj.material.C(:,:,:,:);
            LHS   = zeros(obj.trial.nbasis,obj.test.nbasis);
            for iBasis = 1: obj.trial.nbasis
                uI   = obj.trial.basisFunctions{iBasis};
                defI = uI.computeSymmetricGradient(quad);
                defI = defI.transformInVoigtNotation();
                for jBasis = 1:obj.test.nbasis
                    vJ   = obj.test.basisFunctions{jBasis};
                    defJ = vJ.computeSymmetricGradient(quad);
                    defJ = defJ.transformInVoigtNotation();
                    for kStre = 1:3
                        for lStre = 1:3
                            Ckl = squeeze(Cmat(kStre,lStre,:,:))';
                            ek  = squeeze(defI.fValues(kStre,:,:));
                            el  = squeeze(defJ.fValues(lStre,:,:));
                            kij = ek.*Ckl.*el.*dVolu;
                            LHS(iBasis,jBasis) = LHS(iBasis,jBasis) + sum(kij(:));
                        end
                    end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.test     = cParams.test;
            obj.trial    = cParams.trial;
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.setQuadratureOrder(cParams);
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = 'QUADRATIC';
            end
        end
    end

end