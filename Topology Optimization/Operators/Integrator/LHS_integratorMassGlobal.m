classdef LHS_integratorMassGlobal < handle

    properties (Access = public)

    end

    properties (Access = private)
%         material
        test
        trial
        quadrature
        quadratureOrder
        mesh
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = LHS_integratorMassGlobal(cParams)
            obj.init(cParams)
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            quad  = obj.quadrature;
            posgp = quad.posgp;
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            LHS   = zeros(obj.test.nbasis,obj.trial.nbasis);
            nFlds = obj.trial.ndimf;
            for iBasis = 1: obj.test.nbasis
                uI   = obj.test.basisFunctions{iBasis};
                fI   = uI.evaluate(posgp);
                for jBasis = 1:obj.trial.nbasis
                    vJ   = obj.trial.basisFunctions{jBasis};
                    fJ   = vJ.evaluate(posgp);
                    for iField = 1:nFlds
                        uv = squeeze(fI(iField,:,:).*fJ(iField,:,:));
                        LHS(iBasis,jBasis)  = LHS(iBasis,jBasis) + sum(uv.*dVolu,'all');
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
%             obj.material = cParams.material; %may we need for dynamic
%             problems?
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