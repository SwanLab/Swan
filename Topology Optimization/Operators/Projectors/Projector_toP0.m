classdef Projector_toP0 < Projector

    properties (Access = private)
        quadrature
        M
    end

    methods (Access = public)

        function obj = Projector_toP0(cParams)
            obj.init(cParams);
            obj.computeQuadrature();
            obj.createMassMatrix();
        end

        function xFun = project(obj, x)
            RHS = obj.createRHS(x);
            s.fValues = obj.M\RHS;
            s.mesh    = obj.mesh;
            xFun = P0Function(s);
        end

    end

    methods (Access = private)

        function createMassMatrix(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('CONSTANT');
            dv = obj.mesh.computeDvolume(quad);
            z  = sum(dv(1,:),1);            
            obj.M = spdiags(z', 0, numel(z), numel(z));
        end

        function rhs = createRHS(obj, fun)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV = obj.quadrature.posgp;
            nGaus  = obj.quadrature.ngaus;
            nF     = fun.ndimf;
            nElem  = size(obj.mesh.connec,1);
            rhs = zeros(nElem,nF);
            fGaus = fun.evaluate(xV);

            for iGaus = 1:nGaus
                dVg(:,1) = dV(iGaus,:);
                for iF = 1:nF
                    fGausF = squeeze(fGaus(iF,iGaus,:));
                    Ni = 1;
                    int = Ni*fGausF.*dVg;
                    rhs(:,iF) = rhs(:,iF) + int;
                end
            end
        end

        function computeQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end
    end

end