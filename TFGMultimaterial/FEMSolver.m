classdef FEMSolver < handle

    properties (Access = public)
        charFunc
    end

    properties (Access = private)
        mesh
        pdeCoeff
        tensor
        bc
        E
        gamma
        levelSet
        effectiveTensor
        designVariable
        m % our mesh!!
    end

    methods (Access = public)

        function obj = FEMSolver(cParams)
            obj.init(cParams)
            obj.computeGamma();
            obj.computeCharacteristicFunction();
            obj.computeEffectiveTensor();
        end

        function [U,F] = computeStiffnessMatrixAndForce(obj)
            p = obj.mesh.p;
            t = obj.mesh.t;
            a = obj.pdeCoeff.a;
            f = obj.pdeCoeff.f;
            c = obj.effectiveTensor;
            gamma = obj.gamma;
            tfi = obj.charFunc;
    
          %  tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach
            [K,~,F] = assema(p,t,c,a,f);
            [K,F] = pdeupdate(K,F,obj.bc,obj.mesh);
            U = K \ F;  % solve linear system
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.E(1) = cParams.matProp.matA.young;
            obj.E(2) = cParams.matProp.matB.young;
            obj.E(3) = cParams.matProp.matC.young;
            obj.E(4) = cParams.matProp.matD.young;

            obj.mesh     = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.bc       = cParams.bc;
            obj.levelSet = cParams.psi;
            obj.designVariable = cParams.designVariable;
            obj.m = cParams.m;
        end
        
        function computeGamma(obj)
            obj.gamma = obj.E./obj.E(1); % contrast for each material
        end


        function computeCharacteristicFunction(obj)
            s.p = obj.mesh.p;
            s.t = obj.mesh.t;
            s.pdeCoeff = obj.pdeCoeff; 
            s.psi = obj.levelSet;
            s.designVariable = obj.designVariable;
            s.m = obj.m;

            charfun = CharacteristicFunctionComputer(s); % s'ha de construir la classe - charfunc!!
            [~,tfi] = charfun.computeFiandTfi();
            obj.charFunc = tfi;
        end

        function computeEffectiveTensor(obj)
            tgamma = obj.gamma*obj.charFunc; % for mixed formulation approach
            c0 =  obj.pdeCoeff.tensor(:,1);
            obj.effectiveTensor = c0*tgamma; % effective elasticity tensor
        end

    end

end
