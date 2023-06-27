classdef RHSintegrator_ShapeFunctionWithFunPower < handle

    properties (Access = private)
        mesh
        fun
        exponent
        designVarType
    end

    properties (Access = private)
        interpolationRHS
        quadrature
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeFunctionWithFunPower(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function RHS = compute(obj)
            switch obj.designVarType
                case 'Density'
                    rhs = obj.computeElementalRHS();
                    RHS = obj.assembleIntegrand(rhs);
                case 'LevelSet'

            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh             = cParams.mesh;
            obj.fun              = cParams.fun;
            obj.exponent         = cParams.exponent;
            obj.designVarType    = cParams.designVarType;
            obj.interpolationRHS = Interpolation.create(obj.mesh,'LINEAR');
        end
        
        function createQuadrature(obj)
            quadratureOrder = 'nGaussPointsPerLine';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadratureOrder);
            obj.quadrature = quad;
        end

        function rhs = computeElementalRHS(obj)
            p      = obj.exponent;
            quad   = obj.quadrature;
            dVolu  = obj.mesh.computeDvolume(quad);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nNodes = obj.mesh.nnodeElem;
            nDim   = nNodes*obj.fun.ndimf;
            xg     = quad.posgp;
            obj.interpolationRHS.computeShapeDeriv(xg);
            shapes = obj.interpolationRHS.shape;
            funVals = obj.fun.evaluate(xg);
            rhs = zeros(nDim, 1, nElem);

            for igauss = 1 :nGaus
                for inode= 1:nNodes
                    Ni = shapes(inode,igauss,:);
                    for iunkn= 1:obj.fun.ndimf
                        idof = obj.fun.ndimf*(inode-1)+iunkn;
                        dvol = dVolu(igauss,:);
                        f  = squeeze(funVals(iunkn,igauss,:));
                        rhs(idof, 1, :)= squeeze(rhs(idof,1,:)) ...
                            + ((Ni*f).^p).*dvol';
                    end
                end
            end

        end

        function f = assembleIntegrand(obj,rhsElem)
            integrand = squeeze(rhsElem);
            connec = obj.fun.computeDofConnectivity()';
            ndofs = max(max(connec));
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for iNode = 1:nnode
                int = integrand(iNode,:)';
                con = connec(:,iNode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

    end

end