classdef RHSintegrator_functionWithShapeFunction < handle

    properties (Access = public)
        RHS
    end

    properties (Access = private)
        quadType
        mesh
        fun
        trial
        quadrature
    end

    methods (Access = public)
        function obj = RHSintegrator_functionWithShapeFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeRHS();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadType = cParams.quadType;
            obj.mesh     = cParams.mesh;
            obj.fun      = cParams.fun;
            obj.trial    = cParams.trial;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadType);
            obj.quadrature = q;
        end

        function computeRHS(obj)
            quad = obj.quadrature;
            xV   = quad.posgp;
            dV   = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = obj.trial.computeShapeFunctions(quad);
            nGaus = quad.ngaus;
            nFlds = obj.fun.ndimf;
            nDofElem = size(shapes,1);
            conne = obj.trial.computeDofConnectivity';
            nDofs = max(conne,[],"all");
            fGaus = obj.fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nDofElem
                        dofs = conne(:,inode);
                        Ni = shapes(inode,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            obj.RHS = f;
        end

    end
end