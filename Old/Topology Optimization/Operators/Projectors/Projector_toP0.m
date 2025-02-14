classdef Projector_toP0 < Projector

    properties (Access = private)
        quadrature
        M
        mesh
    end

    methods (Access = public)

        function obj = Projector_toP0(cParams)
          %  obj.init(cParams);
        end

        function xFun = project(obj, x)
            obj.mesh = x.mesh;
            obj.computeQuadrature();  
            obj.createMassMatrix();
            RHS = obj.createRHS(x);
            s.fValues = obj.M\RHS;
            s.mesh    = obj.mesh;
            s.order   = 'P0';
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)

        function createMassMatrix(obj)
            quad = Quadrature.create(obj.mesh,1);
            dv = obj.mesh.computeDvolume(quad);
            a = sum(dv(1,:),1);
            obj.M = spdiags(a',0,length(a),length(a));
         %   obj.M = spdiags(sum(dv(1,:),1),0);
        end

        function rhs = createRHS(obj, fun)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV = obj.quadrature.posgp;
            fGaus = squeezeParticular(fun.evaluate(xV),2);
            nGaus  = obj.quadrature.ngaus;
            nF     = size(fGaus,1);
            nElem  = size(obj.mesh.connec,1);
            rhs = zeros(nElem,nF);

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
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end
    end

end