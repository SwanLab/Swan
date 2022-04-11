classdef NewFilter_P1_Density < handle

    properties (Access = private)
        mesh
        quadratureOrder
        Poper
        x
        x_reg
        geometry
        quadrature
        M
        Kernel
    end


    methods (Access = public)
        
        function obj = NewFilter_P1_Density(cParams)
            obj.init(cParams);
            obj.createMassMatrix(cParams);
            obj.createPoperator(cParams);
            obj.createFilterKernel();
        end

        function preProcess(obj)
            s.mesh            = obj.mesh;
            s.quadratureOrder = obj.quadratureOrder;
            P1proc            = P1preProcessor(s);
            P1proc.preProcess();
            obj.storeParams(P1proc);
        end

        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrateRHS(x0);
            P = obj.Poper.value;
            x_reg = P'*RHS;
        end
        
        function x0 = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                xR = obj.computeP0fromP1(x);
                x0 = zeros(length(xR),obj.quadrature.ngaus);
                for igaus = 1:obj.quadrature.ngaus
                    x0(:,igaus) = xR;
                end
            else
                x0 = obj.x_reg;
            end
            obj.updateStoredValues(x,x0);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createMassMatrix(obj,cParams)
            diffReacProb = obj.createDiffReacProblem(cParams);
            obj.M = diffReacProb.getM();
        end

        function createPoperator(obj,cPar)
            cParams.nelem  = obj.mesh.nelem;
            cParams.nnode  = obj.mesh.nnode;
            cParams.npnod  = obj.mesh.npnod;
            cParams.connec = obj.mesh.connec;
            cParams.diffReactEq = cPar.femSettings;
            obj.Poper = Poperator(cParams);
        end

        function storeParams(obj,P1proc)
            obj.quadrature = P1proc.quadrature;
            obj.geometry   = P1proc.geometry;
        end

        function itHas = xHasChanged(obj,x)
            itHas = ~isequal(x,obj.x);
        end

        function updateStoredValues(obj,x,x0)
            obj.x = x;
            obj.x_reg = x0;
        end

        function x0 = computeP0fromP1(obj,x)
            x0 = obj.Kernel*x;
        end

        function createFilterKernel(obj)
            P = obj.Poper.value;
            obj.Kernel = P*obj.M;
        end

        function intX = integrateRHS(obj,x)
            intX = zeros(obj.mesh.nelem,1);
            ng = size(x,2);
            for igaus = 1:ng
                dvolu = obj.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end

    end

    methods (Access = private, Static)

        function pB = createDiffReacProblem(cParams)
            s = cParams.femSettings;
            s.mesh = cParams.mesh;
            switch s.scale
                case 'MACRO'
                    pB = DiffReactProblem(s);
                case 'MICRO'
                    pB = DiffReactProblemMicro(s);
            end
        end

    end

end