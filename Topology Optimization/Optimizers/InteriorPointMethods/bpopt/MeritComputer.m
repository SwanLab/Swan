classdef MeritComputer < handle
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
        function obj = MeritComputer(cParams)
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
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            u.xL = obj.xL;
            u.xU = obj.xU;
            ph = PhiComputer(u);
            ph.compute();
            obj.phi = ph.phi;
        end

        function computeResidual(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            res = ResidualComputer(u);
            res.compute();
            obj.residual = res.c;
        end

        function computeMeritFunction(obj)
            obj.merit = obj.phi + obj.bp.nu*sum(abs(obj.residual));
        end
    end
end