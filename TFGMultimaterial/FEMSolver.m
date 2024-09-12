classdef FEMSolver < handle

    properties (Access = private)
        mesh
        pdeCoeff
        bc
        effectiveTensor
    end

    methods (Access = public)

        function obj = FEMSolver(cParams)
            obj.init(cParams)
        end

        function [U,F] = computeStiffnessMatrixAndForce(obj)
            p = obj.mesh.p;
            t = obj.mesh.t;
            a = obj.pdeCoeff.a;
            f = obj.pdeCoeff.f;
            c = obj.effectiveTensor;
    
          %  tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach
            [K,~,F] = assema(p,t,c,a,f);
            [K,F] = pdeupdate(K,F,obj.bc,obj.mesh);
            U = K \ F;  % solve linear system
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.effectiveTensor = cParams.effTensor;
            obj.bc       = cParams.bc;
        end

    end

end
