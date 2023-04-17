classdef LHSintegrator_WeakDivergence < handle

    properties (Access = private)
        pressureFun
        velocityFun
        mesh
        quadrature
    end

    methods (Access = public)

        function obj = LHSintegrator_WeakDivergence(cParams)
            %             obj.init(cParams);
            %             obj.createQuadrature();
            %             obj.createInterpolation();
            %             obj.createGeometry();
            obj.initStokesD(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
%             LHS = obj.assembleStokesD(lhs);
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dNdxV = obj.velocityFun.computeCartesianDerivatives(obj.quadrature);
            dvolV = obj.mesh.computeDvolume(obj.quadrature)';
            shpeP = obj.pressureFun.computeShapeFunctions(obj.quadrature);

            nElem = obj.mesh.nelem;
            nDimfV = size(dNdxV,1);
            nNodeV = size(dNdxV,2);
            nNodeP = size(shpeP,1);
            
            nGaus = size(dNdxV,4);

            D = zeros(nDimfV*nNodeV,nNodeP,nElem);
            for igaus = 1:nGaus
                for inode_var = 1:nNodeP
                    for inode_test = 1:nNodeV
                        for idime = 1:nDimfV
                            dof_test = inode_test*nDimfV - nDimfV + idime;
                            v = squeeze(dNdxV(idime,inode_test,:,igaus));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*shpeP(inode_var,igaus)...
                                .*dvolV(:,igaus);
                        end
                    end
                end
            end
            lhs = D;
        end

    end

    methods (Access = private)

        function initStokesD(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.pressureFun = cParams.pressureFun;
            obj.velocityFun = cParams.velocityFun;
%             obj.material = cParams.material;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC'); % ehhh
            obj.quadrature = q;
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assembleFunctions(lhs, obj.velocityFun, obj.pressureFun);
        end

    end

end