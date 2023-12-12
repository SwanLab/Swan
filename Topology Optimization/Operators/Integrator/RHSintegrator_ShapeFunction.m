classdef RHSintegrator_ShapeFunction < handle

    properties (Access = private)
        quadType
        mesh
        quadrature
    end

    methods (Access = public)
        function obj = RHSintegrator_ShapeFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function RHS = compute(obj,fun,test)
            quad = obj.quadrature;
            xV   = quad.posgp;
            dV   = obj.mesh.computeDvolume(quad);
            shapes = test.computeShapeFunctions(quad);
            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nDofElem = size(shapes,1);
            conne = test.computeDofConnectivity';
            nDofs = max(conne,[],"all");
            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for idof = 1:nDofElem
                        dofs = conne(:,idof);
                        Ni = shapes(idof,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
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
