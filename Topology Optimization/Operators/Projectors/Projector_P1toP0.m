classdef Projector_P1toP0 < handle
    
    % Rename: Projector_toP0
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

        function obj = Projector_P1toP0(cParams)
            obj.init(cParams);
            obj.computeQuadrature();
            obj.createMassMatrix();
        end

        function xFun = project(obj, x)
            RHS = obj.createRHS(x);
            s.fValues = obj.M\RHS;
            xFun = P0Function(s);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
        end

        function createMassMatrix(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            dv = obj.mesh.computeDvolume(quad);
            obj.M = diag(dv);
        end

        function rhs = createRHS(obj, fun)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV = obj.quadrature.posgp;
            nGaus  = obj.quadrature.ngaus;
            nF     = size(fun.fValues,1);
            nElem  = size(obj.mesh.connec,1);
            rhs = zeros(nElem,nF);

            % Separate in two loops
            for igaus = 1:nGaus
                xGaus = xV(:,igaus);
                fGaus = fun.evaluate(xGaus);
                dVg(:,1) = dV(igaus,:);
                for iF = 1:nF
                    fGausF = squeeze(fGaus(iF,:,:));
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

