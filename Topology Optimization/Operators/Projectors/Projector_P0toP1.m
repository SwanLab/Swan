classdef Projector_P0toP1 < handle

    properties (Access = private)
        mesh
        connec
        nelem
        nnode
        npnod
        quadOrder
        quadrature
        field
    end
    
    methods (Access = public)

        function obj = Projector_P0toP1(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createField();
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = (LHS\RHS);
            s.type = obj.mesh.type;
            s.connec = obj.mesh.connec;
            s.fNodes = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
            obj.quadOrder = 'LINEAR';
        end
        
        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.field = Field(s);
        end
        
        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fefun)
            xV = obj.quadrature.posgp;
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
            conne = obj.mesh.connec;

            nGaus = obj.quadrature.ngaus;
            nFlds = size(fefun.fElem, 1);
            nElem = size(fefun.fElem, 3);
            nNods = size(shapes,1);
            nNode = size(conne,2);
            nDofs = obj.mesh.nnodes;

            rhs = zeros(nNods,nElem, nFlds);
            f = zeros(nDofs,nFlds);
            for igaus = 1:nGaus
                % fGaus = ...;
                dVg(:,1) = dV(igaus, :);
                for iField = 1:nFlds
                    fG = squeeze(fefun.fElem(iField,:,:));
                    fdVg = fG.*dVg;
                    Ni = shapes(:,:, igaus);
                    rhs(:,:,iField) = rhs(:,:,iField) + bsxfun(@times,Ni,fdVg');
                    for inode = 1:nNode
                        int = rhs(inode,:,iField);
                        con = conne(:,inode);
                        f(:,iField) = f(:,iField) + accumarray(con,int,[nDofs,1],@sum,0);
                    end
                end
            end
            RHS = f;
        end
        
        function shapes = computeShapeFunctions(obj, xG)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(xG);
            shapes = permute(int.shape,[1 3 2]);
        end
        
    end

end

