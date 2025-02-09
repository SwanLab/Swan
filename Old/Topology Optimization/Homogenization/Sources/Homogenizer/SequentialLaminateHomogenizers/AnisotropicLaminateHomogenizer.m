classdef AnisotropicLaminateHomogenizer < handle
    
    properties (Access = private)
        C0
        C1
        theta
        dir
        dim
        qMatrix
        Ltensor
        Stensor
        Rtensor
        ChMean
        ChCorrector
        Ch
        Sh
        P0
        P1
        Itensor
        Ptensor
    end
    
    methods (Access = public)
        
        function obj = AnisotropicLaminateHomogenizer(C0,C1,dir,theta)
            obj.init(C0,C1,dir,theta)
            obj.computeQmatrix();
            obj.computeLtensor();
            obj.computeStensor();
            obj.computeRtensor();
            obj.computeHomogenizedCorrectorTensor();
            obj.computeMeanHomogenizedTensor();
            obj.computeHomogenziedTensor();
            obj.makeHomogenizedTensorVoigt();
            obj.computeAmplificatorTensor()
        end
        
        function C = getHomogenizedTensor(obj)
            C = obj.Ch;
        end
        
        function P = getAmplificatorTensor(obj)
            P = obj.Ptensor;
        end
    end
    
    
    methods (Access = private)
        
        function init(obj,C0,C1,dir,theta)
            obj.C0 = obj.makeItPlaneStress(C0);
            obj.C1 = obj.makeItPlaneStress(C1);
            obj.dir = dir;
            obj.theta = theta;
            obj.dim = obj.C0.getDimension;
        end
        
        
        function C = makeItPlaneStress(obj,C)
            C = Tensor2VoigtConverter.convert(C);
            C = PlaneStressTransformer.transform(C);
            C = Voigt2TensorConverter.convert(C);
        end
        
        function computeQmatrix(obj)
            t  = obj.theta;
            n  = obj.dir.getValue();
            c0 = obj.C0.getValue();
            c1 = obj.C1.getValue();
            d = obj.dim;
            qinv = zeros(d,d);
            for i = 1:d
                for j = 1:d
                    qv = 0;
                    for k = 1:d
                        for l=1:d
                            C0_iklj = c0(i,k,l,j);
                            C1_iklj = c1(i,k,l,j);
                            val = (t*C0_iklj+(1-t)*C1_iklj)*n(k)*n(l);
                            qv = qv + val;
                        end
                    end
                    qinv(i,j) = qv;
                end
            end
            obj.qMatrix = inv(qinv);
        end
        
        function computeLtensor(obj)
            d = obj.dim;
            q = obj.qMatrix;
            I = eye(d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l=1:d
                            obj.Ltensor(i,j,k,l) = 0.5*(q(i,k)*I(j,l) + q(j,k)*I(i,l));
                        end
                    end
                end
            end
            
        end
        
        function computeStensor(obj)
            n  = obj.dir.getValue();
            d  = obj.dim;
            L  = obj.Ltensor;
            S  = zeros(d,d,d,d);
            for m = 1:d
                for i = 1:d
                    for j = 1:d
                        for k = 1:d
                            for l=1:d
                                Lijkl = L(i,j,k,l);
                                S(i,j,m,k) = S(i,j,m,k) + Lijkl*n(l)*n(m);
                            end
                        end
                    end
                end
            end
            obj.Stensor = S;
        end
        
        function computeRtensor(obj)
            S  = obj.Stensor;
            c0 = obj.C0.getValue();
            c1 = obj.C1.getValue();
            R  = obj.computeDoubleProduct(S,c0 - c1);
            obj.Rtensor = R;
        end
        
        function computeHomogenizedCorrectorTensor(obj)
            c0 = obj.C0.getValue();
            c1 = obj.C1.getValue();
            t  = obj.theta;
            R  = obj.Rtensor;
            BAR = obj.computeDoubleProduct(c0-c1,R);
            obj.ChCorrector = t*(t-1)*BAR;
        end
        
        function computeMeanHomogenizedTensor(obj)
            c0 = obj.C0.getValue();
            c1 = obj.C1.getValue();
            t  = obj.theta;
            Chm = (1-t)*c0 + t*c1;
            obj.ChMean = Chm;
        end
        
        function computeHomogenziedTensor(obj)
            obj.Ch = obj.ChMean + obj.ChCorrector;
        end
        
        
        function AB = computeDoubleProduct(obj,A,B)
            d  = obj.dim;
            AB  = zeros(d,d,d,d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l=1:d
                            for m = 1:d
                                for n = 1:d
                                    Aijmn = A(i,j,m,n);
                                    Bmnkl = B(m,n,k,l);
                                    AB(i,j,k,l) = AB(i,j,k,l) + Aijmn*Bmnkl;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function AB = computeTranspouseDoubleProduct(obj,A,B)
            d  = obj.dim;
            AB  = zeros(d,d,d,d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l=1:d
                            for m = 1:d
                                for n = 1:d
                                    Amnij = A(m,n,i,j);
                                    Bmnkl = B(m,n,k,l);
                                    AB(i,j,k,l) = AB(i,j,k,l) + Amnij*Bmnkl;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function makeHomogenizedTensorVoigt(obj)
            Chv = obj.Ch;
            ChT = StiffnessPlaneStressTensor;
            ChT.setValue(Chv);
            obj.Ch = Tensor2VoigtConverter.convert(ChT);
        end
        
        function computeAmplificatorTensor(obj)
            obj.computeComplianceTensor();
            obj.computeIdentityTensor();
            obj.computeP0tensor();
            obj.computeP1tensor();
            P0sq = obj.computeTranspouseDoubleProduct(obj.P0,obj.P0);
            P1sq = obj.computeTranspouseDoubleProduct(obj.P1,obj.P1);
            t = obj.theta;
            p = (1-t)*P0sq + t*P1sq;
            obj.Ptensor = p;
        end
        
        function computeComplianceTensor(obj)
            ch = obj.Ch;
            sh = Inverter.invert(ch);
            obj.Sh = Voigt2TensorConverter.convert(sh);
        end

        function computeP0tensor(obj)
            I = obj.Itensor;
            t = obj.theta;
            R = obj.Rtensor;
            S = obj.Sh.getValue();
            M = obj.computeDoubleProduct(I-(1-t)*R,S);
            c0 = obj.C0.getValue();
            p0 = obj.computeDoubleProduct(c0,M);
            obj.P0 = p0;
        end
        
        function computeP1tensor(obj)
            I = obj.Itensor;
            t = obj.theta;
            R = obj.Rtensor;
            S = obj.Sh.getValue();
            M = obj.computeDoubleProduct(I+(t)*R,S);
            c1 = obj.C1.getValue();
            p1 = obj.computeDoubleProduct(c1,M);
            obj.P1 = p1;
        end
        
        function computeIdentityTensor(obj)
            d = obj.dim;
            I = eye(d);
            Itens = zeros(d,d,d,d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l=1:d
                            Itens(i,j,k,l) = 0.5*(I(i,k)*I(j,l) + I(j,k)*I(i,l));
                        end
                    end
                end
            end
         obj.Itensor = Itens;
        end
    end
    
    
    
    
end

