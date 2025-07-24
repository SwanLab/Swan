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

            extraDim = obj.computeExtraDims(xV);
            nFields  = ndims(fGaus)-extraDim;

            h = zeros(size(fGaus,1:nFields));
            for iField = 1:size(fGaus,1)
                for iGaus = 1:nGaus
                    dVg(:,1) = dV(iGaus, :);
                    if nFields == 1
                        fG = squeeze(fGaus(iField,iGaus,:));
                        int       = fG.*dVg;
                        h(iField) = h(iField) + sum(int);
                    else
                        for jField = 1:size(fGaus,2)
                            fG  = squeeze(fGaus(iField,jField,iGaus,:));
                            int = fG.*dVg;
                            h(iField,jField) = h(iField,jField) + sum(int);
                        end
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
            q = Quadrature.create(obj.mesh,obj.quadType);
            obj.quadrature = q;
        end

        function extraDim = computeExtraDims(obj,xV)
            extraDim = 2;
            if obj.mesh.nelem == 1
                extraDim = extraDim - 1;
                if size(xV,2) == 1
                    extraDim = extraDim -1;
                end
            end
        end

    end

end