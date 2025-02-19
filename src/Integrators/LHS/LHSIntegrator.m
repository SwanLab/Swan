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
 
    end

    methods (Access = protected)

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

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end
    
end