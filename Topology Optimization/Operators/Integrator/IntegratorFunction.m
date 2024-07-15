classdef IntegratorFunction < handle

    properties (Access = private)
        quadType
        mesh
    end

    properties (Access = private)
        quadrature
    end

    methods (Access = public)
        function obj = IntegratorFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function int = compute(obj,f)
            quad      = obj.quadrature;
            xV        = quad.posgp;
            dV        = obj.mesh.computeDvolume(quad);
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
            obj.quadType = cParams.quadType;
            obj.mesh     = cParams.mesh;
        end

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh,obj.quadType);
            obj.quadrature = q;
        end

    end

end