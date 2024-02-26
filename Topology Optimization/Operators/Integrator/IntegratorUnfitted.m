classdef IntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        quadratureInnerMesh
        quadratureInnerCutMesh
    end

    methods (Access = public)
        function obj = IntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createQuadratureInner();
            obj.createQuadratureInnerCut();
        end

        function int = compute(obj,f)
            quad      = obj.quadratureInnerMesh;
            xV        = quad.posgp;
            dV        = obj.unfittedMesh.computeDvolume(quad);
            nGaus     = quad.ngaus;
            fGaus     = f.evaluate(xV);
            nFields   = size(fGaus,1);
            h         = 0;
            for iField = 1:nFields
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG       = squeeze(fGaus(iField,igaus,:));
                    int      = fG.*dVg;
                    h        = h + sum(int);
                end
            end
            int = h;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createQuadratureInner(obj)
            m = obj.unfittedMesh.innerMesh.mesh;
            q = Quadrature.set(m.type);
            q.computeQuadrature(obj.quadType);
            obj.quadratureInnerMesh = q;
        end

        function createQuadratureInnerCut(obj)
            m = obj.unfittedMesh.innerCutMesh.mesh;
            q = Quadrature.set(m.type);
            q.computeQuadrature(obj.quadType);
            obj.quadratureInnerCutMesh = q;
        end
    end

end