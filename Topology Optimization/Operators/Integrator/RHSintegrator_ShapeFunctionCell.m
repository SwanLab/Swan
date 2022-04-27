classdef RHSintegrator_ShapeFunctionCell < handle

    properties (Access = private)
        npnod
        mesh
        globalConnec
        fNodal
        xGauss
        fGauss
        quadOrder
        quadrature
        nunknPerField
    end

    methods (Access = public)

        % Via Integrator_Simple + Integrator
        function obj = RHSintegrator_ShapeFunctionCell(cParams)
            obj.init(cParams);
        end

        function rhs = compute(obj)
            obj.computeQuadrature();
            obj.computeGaussPoints();
            obj.computeFgauss();
            rhs = obj.integrateFgauss();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.fNodal       = cParams.fNodal;
            obj.quadOrder    = cParams.quadOrder;
            obj.globalConnec = cParams.globalConnec;
            obj.nunknPerField = 1;
        end
        
        function computeQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end
        
        function computeGaussPoints(obj)
            q = obj.quadrature;
            xG = repmat(q.posgp,[1,1,obj.mesh.nelem]);
            obj.xGauss = xG;
        end

        function computeFgauss(obj)
            s.fNodes = obj.fNodal;
            s.connec = obj.globalConnec;
            s.type   = obj.mesh.type;
            f = FeFunction(s);
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhs = integrateFgauss(obj)
            rhsCells = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
        function rhsC = computeElementalRHS(obj) % integrate@RHSintegrator
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
%             fdV    = (fG.*dV);
            shapes = obj.computeShapeFunctions();
            nnode  = size(shapes,1);
            nelem  = size(shapes,2);
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                nunkn = obj.nunknPerField;
                for iField = 1:nunkn
                    fdv = fG(igaus,:,iField).*dV(igaus,:);
                    shape = shapes(:, :, igaus);
                    int = int + bsxfun(@times,shape,fdv);
                end
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs = obj.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end
        
        function shapes = computeShapeFunctions(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);
        end

    end

end