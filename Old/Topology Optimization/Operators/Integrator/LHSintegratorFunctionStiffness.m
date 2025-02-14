classdef LHSintegratorFunctionStiffness < handle %LHSintegrator

    properties (Access = private)
        mesh
        test, trial
        quadratureOrder        
        quadrature
        fun
    end

    methods (Access = public)

        function obj = LHSintegratorFunctionStiffness(cParams)
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
            xV = obj.quadrature.posgp;
            dNdxTs = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodETs = size(dNdxTs,2);
            nDofETs = nNodETs*obj.trial.ndimf;
            nNodETr = size(dNdxTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;

            lhs = zeros(nDofETs,nDofETr,nElem);

            fV  = obj.fun.evaluate(xV);
            for igauss = 1 :nGaus
                for inode= 1:nNodETs
                    for jnode= 1:nNodETr
                        for iunkn= 1:obj.test.ndimf
                            for idim = 1:obj.mesh.ndim
                                idof = obj.fun.ndimf*(inode-1)+iunkn;
                                jdof = obj.fun.ndimf*(jnode-1)+iunkn;
                                dvol = dVolu(igauss,:);
                                dNi = dNdxTs(idim,inode,:,igauss);
                                dNj = dNdxTr(idim,jnode,:,igauss);
                                fVG = fV(1,igauss,:);
                                v = squeeze(fVG.*dNi.*dNj);
                                lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                                    + v(:).*dvol';
                            end
                        end
                    end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.fun   = cParams.function;
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

        function Bcomp = createBComputer(obj, fun, dNdx)
            s.fun  = fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputer(s);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end
  
    end

end
