classdef Projector_toP2 < Projector
    
    methods (Access = public)

        function obj = Projector_toP2(cParams)
            obj.init(cParams);
        end

% LHS
%     0.0167   -0.0028   -0.0028   -0.0000   -0.0000   -0.0111
%    -0.0028    0.0167   -0.0028   -0.0000   -0.0111   -0.0000
%    -0.0028   -0.0028    0.0167   -0.0111   -0.0000   -0.0000
%    -0.0000   -0.0000   -0.0111    0.0889    0.0444    0.0444
%    -0.0000   -0.0111   -0.0000    0.0444    0.0889    0.0444
%    -0.0111   -0.0000   -0.0000    0.0444    0.0444    0.0889

% SHAPES
%    -0.0482   -0.0847   -0.0482   -0.0748    0.5176   -0.0748
%    -0.0482   -0.0482   -0.0847   -0.0748   -0.0748    0.5176
%    -0.0847   -0.0482   -0.0482    0.5176   -0.0748   -0.0748
%     0.7955    0.1928    0.1928    0.0335    0.2992    0.2992
%     0.1928    0.7955    0.1928    0.2992    0.0335    0.2992
%     0.1928    0.1928    0.7955    0.2992    0.2992    0.0335
        
        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            xFun = P2Function(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.fun   = P2Function.create(obj.mesh, 1);
            s.quadratureOrder = 'ORDER4';
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            
            f = P2Function.create(obj.mesh, 1);
            shapes = f.computeShapeFunctions(quad);
            conne = f.computeDofConnectivity()';

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nNode = size(conne,2);
            nDofs = max(max(conne));

            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        dofs = conne(:,inode);
                        Ni = shapes(inode,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end