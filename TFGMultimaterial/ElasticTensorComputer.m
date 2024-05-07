classdef ElasticTensorComputer < handle

    properties (Access = public)
        effectiveTensor
        charFunc
    end

    properties (Access = private)
        matProp
        gamma
        mesh
        pdeCoeff
        designVariable
        m
        bc
        E
    end

    methods (Access = public)
        
        function obj = ElasticTensorComputer(cParams)
            obj.init(cParams);
            obj.computeCharacteristicFunction();
            obj.computeGamma();
            obj.computeElasticTensor();
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
            obj.designVariable = cParams.designVariable;
            obj.m = cParams.m;
        end

        function computeCharacteristicFunction(obj)
            s.p = obj.mesh.p;
            s.t = obj.mesh.t;
            s.pdeCoeff = obj.pdeCoeff; 
            s.designVariable = obj.designVariable;
            s.m = obj.m;

            charfun = CharacteristicFunctionComputer(s); % s'ha de construir la classe - charfunc!!
            [~,tfi] = charfun.computeFiandTfi();
            obj.charFunc = tfi;
        end

        function computeGamma(obj)
            obj.gamma = obj.E./obj.E(1); % contrast for each material
        end

        function computeElasticTensor(obj)
            tgamma = obj.gamma*obj.charFunc; % for mixed formulation approach
            c0 =  obj.pdeCoeff.tensor(:,1);
            obj.effectiveTensor = c0*tgamma; % effective elasticity tensor

            % Prova
            cSeba = obj.effectiveTensor;
            lambdaVals = cSeba(6,:);
            muVals     = cSeba(4,:);
            s.order = 'P0';
            s.fValues = lambdaVals;
            s.mesh    = obj.mesh;
            lambdaField = LagrangianFunction(s);
            % ...
            % muField
            % MaterialGiven
            %
        end
    end
end