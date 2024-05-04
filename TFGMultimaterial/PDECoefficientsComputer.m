classdef PDECoefficientsComputer < handle

    properties (Access = public)
        tensor
        a
        f
        b
    end

    properties (Access = private)
      mu
      lambda
      E
    end

    methods (Access = public)

        function obj = PDECoefficientsComputer(cParams)
            obj.init(cParams)
            obj.computeConstitutiveTensor();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mu(1) = cParams.A.shear;
            obj.mu(2) = cParams.B.shear;
            obj.mu(3) = cParams.C.shear;
            obj.mu(4) = cParams.D.shear;

            obj.lambda(1) = cParams.A.lambda;
            obj.lambda(2) = cParams.B.lambda;
            obj.lambda(3) = cParams.C.lambda;
            obj.lambda(4) = cParams.D.lambda;

            obj.E(1) = cParams.A.young;
            obj.E(2) = cParams.B.young;
            obj.E(3) = cParams.C.young;
            obj.E(4) = cParams.D.young;
            
            obj.a = zeros(4,1);
            obj.f = zeros(2,1);
            obj.b = [];
        end

        function computeConstitutiveTensor(obj)
            c=zeros(16,length(obj.E));
            
            c(1,:) = obj.lambda + 2*obj.mu; c(2,:) = 0; c(3,:) = 0; c(4,:) = obj.mu;
            c(5,:) = 0; c(6,:) = obj.lambda; c(7,:) = obj.mu; c(8,:) = 0;
            c(9,:) = c(8,:); c(10,:) = c(7,:); c(11,:) = c(6,:); c(12,:) = c(5,:);
            c(13,:)= c(4,:); c(14,:) = c(3,:); c(15,:) = c(2,:); c(16,:) = c(1,:); 
            % nStre = 3;
            % nGaus = 2;
            % nElem = 25600;
            % 
            % C = zeros(nStre,nStre,nGaus,nElem);
            % C(1,1,:,:)= 2*obj.mu+obj.lambda;
            % C(1,2,:,:)= obj.lambda;
            % C(2,1,:,:)= obj.lambda;
            % C(2,2,:,:)= 2*obj.mu+obj.lambda;
            % C(3,3,:,:)= obj.mu;

            obj.tensor = c;   
        end

    end

end