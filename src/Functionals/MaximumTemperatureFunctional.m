classdef MaximumTemperatureFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        conductivity
        mass
        stateProblem
        quadrature
        saveMax
        saveL2
    end

    methods (Access = public)
        function obj = MaximumTemperatureFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            if size(x,2) == 2
                iter = x{2};
                x=x{1};
            end
            xD      = x.obtainDomainFunction();
            xR      = obj.filterFields(xD);
            xR{1}.setFValues(1 - xR{1}.fValues);          % 1 - FP
            [J, dJ] = obj.computeMaximumTemperatureFunctionAndGradient(xR);
            if iter == 1 || mod(iter,20)== 0
                dJ{1}.print('temp'+string(iter)) 
            end
            if iter ==  600
                max = obj.saveMax; 
                L2 = obj.saveL2;
                save('max.mat','max');
                save('L2.mat','L2');
            end
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.stateProblem = cParams.stateProblem;
            obj.conductivity = cParams.conductivity;
            obj.mass         = cParams.mass;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function u = computeStateVariable(obj, kappa, mass)
            obj.stateProblem.updateConductivity(kappa);
            obj.stateProblem.updateRHSWithMass(mass);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function [J,dJ] = computeMaximumTemperatureFunctionAndGradient(obj, xR)
            kappa  = obj.createDomainFunction(obj.conductivity.fun,xR);           % conductivity on the new domain
            dkappa = obj.createDomainFunction(obj.conductivity.dfun,xR); 
            mass = obj.createDomainFunction(obj.mass.fun,xR);
            u      = obj.computeStateVariable(kappa,mass);                             % solve the PDE
            J = max(u.fValues);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            dJ{1}     = u; % (non diffentiable)
            obj.saveMax(end+1) = J;
            obj.saveL2(end+1) = Norm(u,'L2');
        end

        function f = createDomainFunction(obj,fun,xR)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xR{1},xV);
            s.mesh      = obj.mesh;
            f = DomainFunction(s);
        end

        function fV = createConductivityAsDomainFunction(obj,fun,xR,xV)
            densV = xR.evaluate(xV);
            fV = fun(densV);
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end