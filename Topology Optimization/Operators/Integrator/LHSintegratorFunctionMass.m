classdef LHSintegratorFunctionMass < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
        fun
    end

    methods (Access = public)

        function obj = LHSintegratorFunctionMass(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            quad = obj.quadrature;
            xV   = quad.posgp;
            shapesTest  = obj.test.computeShapeFunctions(xV);
            shapesTrial = obj.trial.computeShapeFunctions(xV);
            dVolu  = obj.mesh.computeDvolume(quad);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);

            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;

            fG = squeeze(obj.fun.evaluate(xV));

            M = zeros(nDofTest, nDofTrial, nElem);
            % lhs = zeros(nDofTest/2, nDofTrial/2, nElem);
            % for igaus = 1:nGaus
            %     dv(1,1,:) = dVolu(igaus,:);
            %     Nv = shapesTest(:,igaus);
            %     Nu = shapesTrial(:,igaus);
            %     NvNu = Nv*Nu';
            %     Aij = bsxfun(@times,NvNu,dv);
            %     lhs = lhs + Aij;
            % end
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iunkn= 1:obj.test.ndimf
                       %     for junkn= 1:obj.trial.ndimf
                                fdv = fG(igauss,:).*dVolu(igauss,:);
                                idof = obj.test.ndimf*(inode-1)+iunkn;
                                jdof = obj.trial.ndimf*(jnode-1)+iunkn;
                                Ni = shapesTest(inode,igauss,:);
                                Nj = shapesTrial(jnode,igauss,:);
                                v = squeeze(Ni.*Nj.*fdv);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:);
                       %     end
                        end
                    end
                end
            end
            lhs = M;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.fun = cParams.function;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.trial.order;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end
end