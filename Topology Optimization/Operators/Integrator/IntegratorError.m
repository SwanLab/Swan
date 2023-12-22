classdef IntegratorError < handle

    properties (Access = private)
        quadType
        mesh
    end

    properties (Access = private)
        quadrature
    end

    methods (Access = public)
        function obj = IntegratorError(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function int = compute(obj,f,g)
            quad      = obj.quadrature;
            xV        = quad.posgp;
            dV        = obj.mesh.computeDvolume(quad);
            nGaus     = quad.ngaus;
            fGaus     = f.evaluate(xV);
            gGaus     = g.evaluate(xV);
            nFields   = size(fGaus,1);
            h         = 0;
            for iField = 1:nFields
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG       = squeeze(fGaus(iField,igaus,:));
                    gG       = squeeze(gGaus(iField,igaus,:));
                    int      = ((fG-gG).^2).*dVg;
                    h        = h + sum(int);
                end
            end
            int = 0.5*h;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadType = cParams.quadType;
            obj.mesh     = cParams.mesh;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadType);
            obj.quadrature = q;
        end

    end

end