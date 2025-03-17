classdef LocalConstraint < handle

    properties (Access = private)
        mesh
        constraint
        unfittedMesh
        epsilon
        target
        filterPer
        value0
    end

    properties (Access = private)
        filterGrad
    end

    methods (Access = public)
        function obj = LocalConstraint(cParams)
            obj.init(cParams);
            obj.filterPer.updateEpsilon(obj.epsilon);
            obj.createFilterAdjoint();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD{1});
            J  = obj.computeFunction(xD{1},xR);
            dJ = obj.computeGradient(xR);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.constraint   = cParams.constraint;
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.epsilon      = cParams.epsilon;
            obj.target       = cParams.target;
            obj.filterPer    = cParams.filterPer;
            obj.value0       = cParams.value0;
        end

        function createFilterAdjoint(obj)
            s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh         = obj.mesh;
            obj.filterGrad = FilterLump(s);
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filterPer.compute(x,2);
        end

        function J = computeFunction(obj,xD,xR)
            difJ = xD.*(1-xR);
            difJ = obj.computeUnfittedFunction(difJ);
            int  = Integrator.compute(difJ,obj.mesh,2);
            J     = 2/(obj.epsilon)*int;
        end

        function fun = computeUnfittedFunction(obj,fun)
            s.uMesh = obj.unfittedMesh;
            s.fun   = fun;
            fun     = UnfittedFunction(s);
        end

        function dJ = computeGradient(obj,xR)
            dj = 2/(obj.epsilon)*(1-2*xR);
            dJ = obj.computeUnfittedFunction(dj);
            dJ = obj.filterGrad.compute(dJ,2);
        end
    end
end