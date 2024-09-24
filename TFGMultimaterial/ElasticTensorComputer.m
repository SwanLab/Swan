classdef ElasticTensorComputer < handle

    properties (Access = public)
        effectiveTensor
        charFunc
        C
        stateProblem
    end

    properties (Access = private)
        matProp
        gamma
        mesh
        pdeCoeff
        designVariable
        m
        bc
        bounCon
        E
        nu
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
            mat = cParams.matProp;
            obj.E(1) = mat.A.young;
            obj.E(2) = mat.B.young;
            obj.E(3) = mat.C.young;
            obj.E(4) = mat.D.young;

            %obj.mesh     = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.bc       = cParams.bc;
            obj.designVariable = cParams.designVariable;
            obj.m = cParams.mesh;
        end

        function computeCharacteristicFunction(obj)
            charfun = obj.designVariable.obtainDomainFunction();
            [~,tfi] = charfun.computeAtNodesAndElements();
            obj.charFunc = tfi;
        end

        function computeGamma(obj)
            obj.gamma = obj.E./obj.E(1); % contrast for each material
        end

        function computeElasticTensor(obj)
            chi = obj.charFunc.fValues';
            tgamma = obj.gamma*chi; % for mixed formulation approach
            c0 =  obj.pdeCoeff.tensor(:,1);
            obj.effectiveTensor = c0*tgamma; % effective elasticity tensor
            
            % Swan Elastic Tensor
            cSeba = obj.effectiveTensor;
            lambdaVals = cSeba(6,:);
            muVals     = cSeba(4,:);
            
            s.order = 'P0';
            s.fValues = lambdaVals';
            s.mesh    = obj.m;
            lambdaField = LagrangianFunction(s);

            s.order = 'P0';
            s.fValues = muVals';
            s.mesh = obj.m;
            muField = LagrangianFunction(s);

            N = obj.m.ndim;
            kField = lambdaField + 2.*muField/N;

            % Material Given
            s.type  = 'ISOTROPIC';
            s.ptype = 'ELASTIC';
            s.ndim  = N;
            s.shear = muField;
            s.bulk  = kField;
            obj.C    = Material.create(s);
        end

    end
end