classdef PerimeterNormPFunctional < handle

    properties (Access = private)
        quadrature
        Pp
        totalPerimeter
        filter
        epsilon
    end

    properties (Access = private)
        mesh
        perimeterTarget
        p
    end

    methods (Access = public)
        function obj = PerimeterNormPFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalPerimeter();
            obj.createFilter();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)   % dJ: integrand field
            Le    = obj.computeFilteredVariable(x); % rho
            J     = obj.computeFunction(x,Le);
            dJ{1} = obj.computeGradient(x,Le);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh            = cParams.mesh;
            obj.perimeterTarget = cParams.perimeterTarget;
            obj.p               = cParams.p;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,3);
            obj.quadrature = quad;
        end

        function createTotalPerimeter(obj) % revisar
            a = obj.mesh.coord(end,1);
            b = obj.mesh.coord(end,2);
            obj.totalPerimeter = 2*a + 2*b;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.epsilon  = 6*obj.mesh.computeMeanCellSize();
            f.updateEpsilon(obj.epsilon);
            obj.filter = f;
        end

        function xD = computeFilteredVariable(obj,x)
            xD = x.obtainDomainFunction();
            xD = obj.filter.compute(xD{1},3);
        end

        function J = computeFunction(obj,x,Le)
            xP     = ((1-Le)*x)^obj.p;
            PerP   = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            obj.Pp = 1/(2*obj.epsilon)*PerP^(1/obj.p);
            J      = ((1/obj.perimeterTarget)*((1/obj.totalPerimeter)^(1/obj.p))*obj.Pp) - 1;
        end

        function dJ = computeGradient(obj,x,Le)
            num = (obj.Pp^(1-obj.p))*(((1-Le)*x)^(obj.p-1))*(1-2*Le);
            den = 2*obj.epsilon*obj.perimeterTarget*(obj.totalPerimeter)^(1/obj.p);
            dJ  = num/den;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter p-norm';
        end
    end
end