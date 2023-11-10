classdef IntegratorEnergy < handle

    properties (Access = private)
        quadType
        mesh
    end

    properties (Access = private)
        quadrature
    end

    methods (Access = public)
        function obj = IntegratorEnergy(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function int = compute(obj,e,C)
            quad      = obj.quadrature;
            xV        = quad.posgp;
            dV        = obj.mesh.computeDvolume(quad);
            nGaus     = quad.ngaus;
            eGaus     = e.evaluate(xV);
            nFields   = size(eGaus,1);
            h         = 0;
            for iField = 1:nFields
                for jField = 1:nFields
                    for igaus = 1:nGaus
                        dVg(:,1) = dV(igaus, :);
                        energy = eGaus(iField,igaus,:).*C(iField,jField,:,igaus).*eGaus(jField,igaus,:);
                        energyG       = squeeze(energy);

                        int = energyG.*dVg;
                        h   = h + sum(int);
                    end
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
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadType);
            obj.quadrature = q;
        end

    end

end