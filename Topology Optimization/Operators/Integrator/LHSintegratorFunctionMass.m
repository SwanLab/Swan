classdef LHSintegratorFunctionMass < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegrator_Mass(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS(fun);
            LHS = obj.assembleMatrix(fun,lhs);
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
                obj.quadratureOrder = obj.trial.order;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function lhs = computeElementalLHS(obj,fun)
            quad = obj.quadrature;
            shapesTest  = obj.test.computeShapeFunctions(quad);
            shapesTrial = obj.trial.computeShapeFunctions(quad);
            dVolu  = obj.mesh.computeDvolume(quad);
            fG = squeeze(fun.evaluate(quad.posgp));
            
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;

            M = zeros(nDofTest, nDofTrial, nElem);
%             for igaus = 1:nGaus
%                 dv(1,1,:) = dVolu(igaus,:);
%                 Nv = shapesTest(:,igaus);
%                 Nu = shapesTrial(:,igaus);
%                 NvNu = Nv*Nu';
%                 Aij = bsxfun(@times,NvNu,dv);
%                 lhs = lhs + Aij;
%             end
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
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*fdv;
                       %     end
                        end
                    end
                end
            end
            lhs = M;

        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assembleFunctions(lhs, obj.test, obj.trial);
        end

    end
end