classdef LHSintegrator_MassTestTrial < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegrator_MassTestTrial(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
%                 obj.quadratureOrder = obj.test.order;
                obj.quadratureOrder = 'LINEAR';
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function lhs = computeElementalLHS(obj)
            quad = obj.quadrature;
            shapesTest  = obj.test.computeShapeFunctions(quad);
            shapesTrial = obj.trial.computeShapeFunctions(quad);
            dVolu  = obj.mesh.computeDvolume(quad);

            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nDimf  = obj.test.ndimf;
%             nDofs  = numel(obj.fun.fValues);
            nNodE  = size(shapesTest,1);
            nDofE  = nNodE*nDimf;

            lhs = zeros(size(shapesTrial,1), size(shapesTest,1), nElem);
            for igaus = 1:nGaus
                dv(1,1,:) = dVolu(igaus,:);
                Nv = shapesTest(:,igaus);
                Nu = shapesTrial(:,igaus);
                NvNu = Nv*Nu';
                Aij = bsxfun(@times,NvNu,dv);
                lhs = lhs + Aij;
            end

        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assembleFunctions(lhs, obj.test, obj.trial);
        end

    end

end