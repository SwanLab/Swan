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
            LHS   = zeros(obj.trial.nbasis,obj.test.nbasis);
            LHS2   = zeros(obj.trial.nbasis,obj.test.nbasis);
            nFlds = obj.test.ndimf;
            basisTest  = obj.test.evaluateBasisFunctions(posgp);
            basisTrial = basisTest;
            for iBasis = 1: obj.trial.nbasis
                uI   = obj.trial.basisFunctions{iBasis};
                fI   = uI.evaluate(posgp);
                for jBasis = 1:obj.test.nbasis
                    vJ   = obj.test.basisFunctions{jBasis};
                    fJ   = vJ.evaluate(posgp);
                    for iField = nFlds
                        uv = squeeze(fI(iField,:,:).*fJ(iField,:,:));
                        uv2= squeeze(basisTest{iBasis}(iField,:,:).*basisTrial{jBasis}(iField,:,:));
                        
                        LHS(iBasis,jBasis)  = LHS(iBasis,jBasis) + sum(uv.*dVolu,'all');
                        LHS2(iBasis,jBasis)  = LHS2(iBasis,jBasis) + sum(uv2.*dVolu,'all');
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