classdef PDECoefficientsComputer < handle

    properties (Access = public)
        tensor
        a
        f
        b
        C
    end

    properties (Access = private)
      mu
      lambda
      E
      mesh
    end

    methods (Access = public)

        function obj = PDECoefficientsComputer(cParams)
            obj.init(cParams)
            obj.computeConstitutiveTensor();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            mat = cParams.mat;
            obj.mu(1) = mat.A.shear;
            obj.mu(2) = mat.B.shear;
            obj.mu(3) = mat.C.shear;
            obj.mu(4) = mat.D.shear;

            obj.lambda(1) = mat.A.lambda;
            obj.lambda(2) = mat.B.lambda;
            obj.lambda(3) = mat.C.lambda;
            obj.lambda(4) = mat.D.lambda;

            obj.E(1) = mat.A.young;
            obj.E(2) = mat.B.young;
            obj.E(3) = mat.C.young;
            obj.E(4) = mat.D.young;
            
            obj.a = zeros(4,1);
            obj.f = zeros(2,1);
            obj.b = [];
            obj.mesh = cParams.m;
        end

        function computeConstitutiveTensor(obj)
            c=zeros(16,length(obj.E));
            
            c(1,:) = obj.lambda + 2*obj.mu; c(2,:) = 0; c(3,:) = 0; c(4,:) = obj.mu;
            c(5,:) = 0; c(6,:) = obj.lambda; c(7,:) = obj.mu; c(8,:) = 0;
            c(9,:) = c(8,:); c(10,:) = c(7,:); c(11,:) = c(6,:); c(12,:) = c(5,:);
            c(13,:)= c(4,:); c(14,:) = c(3,:); c(15,:) = c(2,:); c(16,:) = c(1,:); 
            
            obj.tensor = c;

            % Mat 1
            s.lambdaField = AnalyticalFunction.create(@(x) obj.lambda(1)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.muField     = AnalyticalFunction.create(@(x) obj.mu(1)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.type        = 'Given';
            tensor1       = Material.create(s);
            
            quad          = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            quadrature    = quad;
            xV            = quadrature.posgp;
            obj.C{1}    = tensor1.evaluate(xV);
            
            % Mat 2
            s.lambdaField = AnalyticalFunction.create(@(x) obj.lambda(2)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.muField     = AnalyticalFunction.create(@(x) obj.mu(2)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.type        = 'Given';
            tensor2       = Material.create(s);
            
            quad          = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            quadrature    = quad;
            xV            = quadrature.posgp;
            obj.C{2}    = tensor2.evaluate(xV);

            % Mat 3
            s.lambdaField = AnalyticalFunction.create(@(x) obj.lambda(3)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.muField     = AnalyticalFunction.create(@(x) obj.mu(3)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.type        = 'Given';
            tensor3       = Material.create(s);
            
            quad          = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            quadrature    = quad;
            xV            = quadrature.posgp;
            obj.C{3}    = tensor3.evaluate(xV);

            % Mat 4
            s.lambdaField = AnalyticalFunction.create(@(x) obj.lambda(4)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.muField     = AnalyticalFunction.create(@(x) obj.mu(4)*ones(size(x(1,:,:))), 2, obj.mesh);
            s.type        = 'Given';
            tensor4       = Material.create(s);
            
            quad          = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            quadrature    = quad;
            xV            = quadrature.posgp;
            obj.C{4}    = tensor4.evaluate(xV);

            
        end

    end

end