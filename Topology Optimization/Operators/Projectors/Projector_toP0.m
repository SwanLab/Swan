classdef Projector_toP0 < Projector
    
    properties (Access = public)
    end

    properties (Access = private)
        mesh
        connec
    end

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
            s.connec  = obj.connec;
            s.type    = obj.mesh.type;
            xFun = P0Function(s);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
        end

        function createMassMatrix(obj)
            % No connectivities (1,2,3,4,..)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            dv = obj.mesh.computeDvolume(quad);
            obj.M = diag(dv(1,:));
        end

        function rhs = createRHS(obj, fun)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV = obj.quadrature.posgp;
            nGaus  = obj.quadrature.ngaus;
            nF     = fun.ndimf;
            nElem  = size(obj.mesh.connec,1);
            rhs = zeros(nElem,nF);
            fGaus = fun.evaluate(xV);

            % Separate in two loops
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

