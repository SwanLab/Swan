classdef Filter_P1_Density < handle

    properties (Access = private)
        Poper
        x
        x_reg
        M
        Kernel
        field
        fieldM
    end

    properties (Access = private)
        mesh
        quadratureOrder
    end

    methods (Access = public)

        function obj = Filter_P1_Density(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            obj.createPoperator(cParams);
            obj.createFilterKernel();
        end

        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrateRHS(x0);
            P = obj.Poper.value;
            x_reg = P'*RHS;
        end
        
        function x0 = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                xR = obj.computeP0fromP1(x);
                x0 = zeros(length(xR),obj.field.quadrature.ngaus);
                for igaus = 1:obj.field.quadrature.ngaus
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
            obj.createField();
            obj.createFieldMass();
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'LINEAR';
            obj.field = Field(s);
        end

        function createFieldMass(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.fieldM = Field(s);
        end

        function createMassMatrix(obj)
            s.dim          = obj.fieldM.dim;
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.field        = obj.fieldM;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function createPoperator(obj,cPar)
            cParams.nelem  = obj.mesh.nelem;
            cParams.nnode  = obj.mesh.nnodeElem;
            cParams.npnod  = obj.mesh.nnodes;
            cParams.connec = obj.mesh.connec;
            cParams.diffReactEq = cPar.femSettings;
            obj.Poper = Poperator(cParams);
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
                dvolu = obj.field.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end

    end

end