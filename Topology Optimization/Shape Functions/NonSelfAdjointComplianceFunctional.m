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
        material
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
            [C,dC]     = obj.computeTensorFunctionAndGradient(xR);
            uS         = obj.computeStateVariable(C);
            uA         = obj.computeAdjointVariable(C);
            J          = obj.computeFunctionValue(C,uS,uA);
            dJ         = obj.computeGradient(dC,uS,uA);
            dJ         = obj.filter.compute(dJ,'LINEAR');
            dJ.fValues = dJ.fValues/obj.value0;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.material     = cParams.material;
            obj.stateProblem = cParams.stateProblem;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');
            obj.quadrature = quad;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,'LINEAR');
        end

        function [C,dC] = computeTensorFunctionAndGradient(obj,xR)
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
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
            s.boundaryConditions.pointloadFun   = bcAdj.pointloadFun;
            s.boundaryConditions.pointload_dofs = bcAdj.pointload_dofs;
            s.boundaryConditions.pointload_vals = bcAdj.pointload_vals;
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