classdef primaldualdeblur < handle
    
    
    properties (Access = protected)
        u0
        M
        N
        MN
        D
        h
        gD
        L2
        tau
        sig
        lambda
        En
    end
    
    
    methods (Access = public)
        
        function obj = primaldualdeblur()
            obj.readImage();
            obj.createDmatrix();
            obj.createH();
            obj.createGDvector();
            obj.plotFirstImage();
            obj.computeParameters()
            obj.compute();
        end
    end
    
    methods (Access = private)
        
        function readImage(obj)
            u0 = double(imread('data/einstein.png'));
            obj.u0 = mean(u0,3);            
            [M N] = size(u0);
            MN = M*N;
            obj.M = M;
            obj.N = N;
            obj.MN = MN;
        end
        
        function createDmatrix(obj)
            M = obj.M;
            N = obj.N;
            MN = obj.MN;
            I=reshape([1:M*N],M,N);
            east=[I(:,2:end), I(:,end)];
            north=[I(2:end,:); I(end,:)];
            D1 = sparse(I,east,1,MN,MN) -speye(MN,MN);
            D2 = sparse(I,north,1,MN,MN) -speye(MN,MN);
            obj.D = [D1 ; D2];
        end
        
        function createH(obj)
            h0 = [ 1 2 1 ];
            h = h0;
            for i=1:1 %% level of blur
                h = conv(h0,h);
            end
            h = h'*h;
            h = h/sum(sum(h));
            obj.h = h;
        end
        
        function createGDvector(obj)
            h = obj.h;
            u0 = obj.u0;
            gD = obj.gD;
            stddev=2;
            gD = filter2(h,u0,'valid');
            [Ms Ns]=size(gD);
            obj.gD = gD + stddev*randn(Ms,Ns);
        end
        
        function plotFirstImage(obj)
            figure(1);
            imagesc(obj.u0); colormap(gray); %drawnow();
            figure(2);
            imagesc(obj.gD); colormap(gray); drawnow();
        end
        
        function computeParameters(obj)
            obj.L2 = 8;
            obj.tau = 1; %% max value for tau = 1;
            obj.lambda = .1; %lambda = .2;
            obj.sig = 1/(obj.tau*obj.L2);
        end
        
        function compute(obj)
            D = obj.D;
            tau = obj.tau;
            sig = obj.sig;
            lambda = obj.lambda;
            
            u=zeros(size(obj.u0));            
            p = zeros(2*obj.MN,1);
            %% primal-dual
            for i=1:1000
                um = u;
                grad = obj.computeGradient(p,u);
                u = u-tau*(grad);
                um = 2*u-um;
                p = p + sig*D*um(:);    
                
                
              %  p = obj.computePdirection(u,p);
              %  u = obj.computeU(u,p);
                
                
                
                p = obj.projectInTheBall(p);
                obj.computeEnergy(i,p,u);
                obj.printImage(i,u);
            end

        end
        
        
        function p = computePdirection(obj,u,p)
            D = obj.D;
            g = obj.computeG(u);
            tau = obj.tau;
            p = p - tau*D*(D'*p - g);
        end
        
        function u = computeU(u,p)
            g = obj.computeG(u);
            D = obj.D;
            u = g(:) - D'*p;
        end
        
        function normP = computeNormP(obj,p)
             MN = obj.MN;
             normP = hypot(p(1:MN),p(MN+1:end));            
        end
        
        function computeEnergy(obj,i,p,u)
            D = obj.D;
            g = obj.computeG(u);
            obj.En(i) = sum ((D'*p).^2)/2 - sum((D'*p).*g(:));
        end
        
        function p = projectInTheBall(obj,p)
            lambda = obj.lambda;
            normP = obj.computeNormP(p);
            no= max(1,normP/lambda);
            p = p./[no;no];
        end
        
        function grad = computeGradient(obj,p,u)
            D = obj.D;
            M = obj.M;
            N = obj.N;
            g = obj.computeG(u);
            grad = reshape(D'*p,M,N) - g;            
        end
        
        function g = computeG(obj,u)
            h = obj.h;
            gD = obj.gD;
            g = conv2(h,filter2(h,u,'valid')-gD,'full');
        end
        
        function printImage(obj,i,u)
            if mod(i,20)==0
                i
                figure(3);
                imagesc(u);
                colormap(gray); drawnow();
            end
        end
        
    end
    
    
    
    
end
