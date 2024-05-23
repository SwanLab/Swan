classdef LocalPerimeterFunctional < handle

    properties (Access = private)
        uMesh
        filter
        epsilon
        value0
        subDomain
        old2new
    end

    methods (Access = public)
        function obj = LocalPerimeterFunctional(cParams)
            obj.init(cParams);
            obj.filter.updateEpsilon(obj.epsilon);
            obj.createSubDomainMesh(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            xD = obj.uMesh.obtainFunctionAtUnfittedMesh(xD);
            xRei = obj.uMesh.obtainFunctionAtUnfittedMesh(xR);
            xD = xD.innerMeshFunction;
            xR = xRei.innerMeshFunction;
            J  = obj.computeFunction(xD,xR);
            dJ = obj.computeGradient(xRei);
            J  = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);
        end

        function updateEpsilon(obj,epsilon)
            obj.filter.updateEpsilon(epsilon);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.filter    = cParams.filter;
            obj.epsilon   = cParams.epsilon;
            obj.value0    = cParams.value0;
        end

        function createSubDomainMesh(obj,cParams)
            fH                 = cParams.subDomainHandle;
            sLS.type           = 'Given';
            sLS.fHandle        = fH;
            g                  = GeometricalFunction(sLS);
            lsFun              = g.computeLevelSetFunction(cParams.mesh);
            levelSet           = lsFun.fValues;
            sUm.backgroundMesh = cParams.mesh;
            sUm.boundaryMesh   = cParams.mesh.createBoundaryMesh();
            obj.uMesh          = UnfittedMesh(sUm);
            obj.uMesh.compute(levelSet);
            [obj.subDomain,obj.old2new] = obj.uMesh.createInnerMesh();
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
        end

        function J = computeFunction(obj,xD,xR)
            f   = xD.*(1-xR);
            int = Integrator.compute(f,obj.uMesh.innerMesh.mesh,2);
            J   = 2/(obj.epsilon)*int;
        end

        function dJ = computeGradient(obj,xRei)
%             dj         = 2/(obj.epsilon)*(1-2*xR.fValues);
%             dJ         = LagrangianFunction.create(obj.subDomain,1,'P1');
%             dJ.fValues = dj;

            M = obj.createMassMatrix();

            s.mesh = obj.uMesh.innerMesh.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 2;
            int        = RHSintegrator.create(s);
            test       = LagrangianFunction.create(s.mesh,1,'P1');
            rhsI        = int.compute(xRei.innerMeshFunction,test); % dJ

            s.mesh = obj.uMesh.innerCutMesh.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 2;
            int        = RHSintegrator.create(s);
            test       = LagrangianFunction.create(s.mesh,1,'P1');
            rhsIC        = int.compute(xRei.innerCutMeshFunction,test); % dJ

            rhs = [rhsI;rhsIC];
            newRHS = M'*rhs;
            isZero = newRHS==0;
            LHS = M'*M;
            LHSfree = LHS(~isZero,~isZero);
            fValuesfree=newRHS(~isZero)./sum(LHSfree,2);
            newfValues = zeros(size(newRHS));
            newfValues(~isZero) = fValuesfree;
            xR2=LagrangianFunction.create(obj.uMesh.backgroundMesh,1,'P1');
            xR2.fValues = newfValues;
        end

        function M = createMassMatrix(obj)
            s.uMesh       = obj.uMesh;
            s.testOrder   = 'P1';
            s.trialOrder  = 'P1';
            s.quadratureOrder = 2;
            s.type  = 'MassMatrixDiffMeshes';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Loc. perimeter';
        end
    end
end