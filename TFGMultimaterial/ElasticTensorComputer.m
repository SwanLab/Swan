classdef ElasticTensorComputer < handle

    properties (Access = public)
        effectiveTensor
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
        nu1
    end

    methods (Access = public)
        
        function obj = ElasticTensorComputer(cParams)
            obj.init(cParams);
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
            obj.nu1  = mat.A.nu;

            %obj.mesh     = cParams.mesh;
            obj.pdeCoeff = cParams.pdeCoeff;
            obj.bc       = cParams.bc;
            obj.designVariable = cParams.designVariable;
            obj.m = cParams.mesh;
        end

        function computeGamma(obj)
            obj.gamma = obj.E./obj.E(1); % contrast for each material
        end

        function computeElasticTensor(obj)
            x = obj.designVariable;
            c0 =  obj.pdeCoeff.tensor(:,1);
            
            sC.E = obj.E;
            sC.C0 = c0;
            sC.nu1 = obj.nu1;
            intMat = MultiMaterialInterpolation(sC);
            [mu,kappa] = intMat.computeConsitutiveTensor(x);
            
            N = obj.m.ndim;
            % Material Given
            s.type  = 'ISOTROPIC';
            s.ptype = 'ELASTIC';
            s.ndim  = N;
            s.shear = mu;
            s.bulk  = kappa;
            obj.C    = Material.create(s);
        end

    end
end