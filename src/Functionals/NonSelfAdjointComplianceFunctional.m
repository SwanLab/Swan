classdef NonSelfAdjointComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        quadrature
        adjointProblem
    end

    properties (Access = private)
        mesh
        filter
        C
        dC
        stateProblem
    end

    methods (Access = public)
        function obj = NonSelfAdjointComplianceFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createAdjointProblem(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD         = x.obtainDomainFunction();
            xR         = obj.filterDesignVariable(xD);
            [CxR,dCxR] = obj.computeTensorFunctionAndGradient(xR);
            uS         = obj.computeStateVariable(CxR);
            uA         = obj.computeAdjointVariable(CxR);
            J          = obj.computeFunctionValue(CxR,uS,uA);
            dJ         = obj.computeGradient(dCxR{1},uS,uA);
            dJ         = {obj.filter.compute(dJ,2)};
            dJVal      = dJ{1}.fValues/obj.value0;
            dJ{1}.setFValues(dJVal);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.C            = cParams.C;
            obj.dC           = cParams.dC;
            obj.stateProblem = cParams.stateProblem;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, 2);
            obj.quadrature = quad;
        end

        function xR = filterDesignVariable(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [CxR,dCxR] = computeTensorFunctionAndGradient(obj,xR)
            CxR     = obj.C(xR{1});
            dCxR{1} = obj.dC(xR{1});
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function u = computeAdjointVariable(obj,C)
            obj.adjointProblem.updateMaterial(C);
            obj.adjointProblem.solve();
            u = obj.adjointProblem.uFun;
        end

        function J = computeFunctionValue(obj,C,uS,uA)
            stateStrain   = SymGrad(uS);
            adjointStrain = SymGrad(uA);
            stressA       = DDP(C,adjointStrain);
            dCompliance   = DDP(stateStrain,stressA);
            J             = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;
        end

        function createAdjointProblem(obj,cParams)
            file                 = cParams.filename;
            [fAdj, fAdj2]        = Preprocess.getBC_adjoint(file, obj.mesh);
            a.fileName           = file;
            s                    = FemDataContainer(a);
            s.bc.pointload       = fAdj;
            s.newBC.pointloadFun = fAdj2;
            bcAdj = obj.getAdjointBoundaryConditions(fAdj2);
            s.boundaryConditions.tractionFun   = bcAdj.tractionFun;
            obj.adjointProblem = PhysicalProblem.create(s);
        end
        
        function bc = getAdjointBoundaryConditions(obj, fAdj2)
            a.mesh         = obj.mesh;
            a.pointloadFun = fAdj2;
            a.dirichletFun = [];
            a.periodicFun  = [];
            bc             = BoundaryConditions(a);
        end
    end

    methods (Static, Access = private)
        function dJ = computeGradient(dC,uS,uA)
            stateStrain   = SymGrad(uS);
            adjointStrain = SymGrad(uA);
            dStressA      = DDP(dC,adjointStrain);
            dJ            = -DDP(stateStrain,dStressA);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'NonSelfAdjCompliance';
        end
    end
end