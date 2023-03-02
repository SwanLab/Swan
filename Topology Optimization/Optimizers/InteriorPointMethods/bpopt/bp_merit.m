classdef bp_merit < handle
    properties (Access = public)
        merit
    end
    properties (Access = private)
        phi
        residual
        x 
        s 
        xL 
        xU 
        bL 
        bU 
        bp 
    end

    methods (Access = public)
        function obj = bp_merit(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computePhi();
            obj.computeResidual();
            obj.computeMeritFunction();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.xL = cParams.xL;
            obj.xU = cParams.xU;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
            obj.bp = cParams.bp;
        end

        function computePhi(obj)
            ph = bp_phi(obj);
            ph.compute();
            obj.phi = ph.phi;
        end

        function computeResidual(obj)
            res = bp_res(obj);
            res.compute();
            obj.residual = res.c;
        end

        function computeMeritFunction(obj)
            obj.merit = obj.phi + obj.bp.nu*sum(abs(obj.residual));
        end
    end
end

function [me] = bp_merit(bp,x,xL,xU,s,bL,bU)
    ph = bp_phi(bp,x,xL,xU,s,bL,bU);
    r = bp_res(bp,x,s,bL,bU);
    me = ph + bp.nu*sum(abs(r));