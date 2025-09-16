classdef LHSIntegrator < handle

    properties (Access = protected)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSIntegratorFactory();
            obj = f.create(s);
        end

    end

    methods (Access = public)
        
        function obj = LHSIntegrator(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj,f,test,trial)
            lhs = obj.computeElementalLHS(f);
            LHS = obj.assembleMatrix(lhs,test,trial);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj,f)
            nElem  = obj.mesh.nelem;
            lhs    = zeros(size(f,1),size(f,2),nElem);
            J = Jacobian(obj.mesh);
            detJ = Det(J);

            xV = obj.quadrature.posgp;
            w  = obj.quadrature.weigp;
            for i = 1:size(f,1)
                for j = 1:size(f,2)
                    int = (f{i,j}.*detJ)*w';
                    lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
                end
            end
        end

        function init(obj, cParams)
            %obj.test  = cParams.test;
            %obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                % warning('Assuming quadrature order')
                quadOrderTe = obj.test.getOrderNum();
                quadOrderTr = obj.trial.getOrderNum();
                obj.quadratureOrder = quadOrderTe + quadOrderTr;
            end
        end
        
        function createQuadrature(obj)
%             quadOrder = obj.fun.getOrderNum();
            quad = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function LHS = assembleMatrix(obj, lhs,test,trial)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, test, trial);
        end

    end
    
end