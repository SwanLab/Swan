classdef ForwardBackwardExample < handle
    
    
    properties (Access = private)
        u0
        rows
        columns
        rowsColumns
        Differential
        
        betaAll
        
        uStar
        L2
        tau
        sig
        lambda
        En
        EnA
        tOld
        tOldOld
        pOldOld
        maxIter
        u
        pDual
        gap
        gapA
        gNorm
    end
    
    
    methods (Access = public)
        
        function obj = ForwardBackwardExample()
            obj.readImage();
            obj.createDmatrix();
            obj.computeParameters();
            obj.computeG();
            obj.computeGnorm();
            obj.plotInitialImage();
            obj.computeAccelerated();
            obj.compute();
            obj.plotEnergy();
            obj.plotDualGap();
        end
        
    end
    
    
    methods (Access = private)
        

        function readImage(obj)
            uInit = double(imread('data/einstein.png'));
            [m,n] = size(uInit);
            obj.u0 = uInit(:); 
            obj.rows = m;
            obj.columns = n;
            obj.rowsColumns = m*n;
        end
        
        function createDmatrix(obj)
            M = obj.rows;
            N = obj.columns;
            MN = obj.rowsColumns;
            I=reshape(1:M*N,M,N);
            east=[I(:,2:end), I(:,end)];
            north=[I(2:end,:); I(end,:)];
            D1 = sparse(I,east,1,MN,MN) -speye(MN,MN);
            D2 = sparse(I,north,1,MN,MN) -speye(MN,MN);
            obj.Differential = [D1 ; D2];
        end
        
        function computeParameters(obj)
            obj.L2 = 8;
            obj.tau = 1/(obj.L2); 
            obj.lambda = 10;
            obj.sig = 100/(obj.L2);
            obj.maxIter = 10;
        end
        
        function computeG(obj)
            u0v = obj.u0;
            sigma = obj.sig;
            obj.uStar= u0v+sigma*randn(size(u0v));
        end
        
        function computeGnorm(obj)
            g = obj.uStar;
            obj.gNorm = 0.5*g'*g;
        end
        
        function compute(obj)
            obj.pDual = zeros(2*obj.rowsColumns,1);
            for i=1:obj.maxIter
                obj.computeGradientStep();
                obj.projectInTheBall();
                obj.computeU();                
                obj.En(i) = obj.computeEnergy();
                obj.gap(i) = obj.computeDualGap();
                obj.printIterationImage(i);
            end

        end
        
        function computeAccelerated(obj)
            obj.pOldOld = zeros(2*obj.rowsColumns,1);
            obj.pDual   = zeros(2*obj.rowsColumns,1);
            for i=1:obj.maxIter
                pOld = obj.pDual;
                obj.computeBetaStep(i);
                obj.computeGradientStep();
                obj.projectInTheBall();
                obj.pOldOld = pOld;
                obj.computeU();                
                obj.EnA(i)  = obj.computeEnergy();
                obj.gapA(i) = obj.computeDualGap();
                obj.printIterationImage(i);
            end
        end
        
        function computeBetaStep(obj,i)
            pOld = obj.pOldOld;
            p    = obj.pDual;
            beta = obj.computeBeta(i);
            obj.betaAll(i) = beta;
            p = p + beta*(p - pOld);
            obj.pDual = p;
        end
        
        function beta = computeBeta(obj,i)
            if i < 3
                obj.tOld = 1;
                obj.tOldOld = 1;
                beta = 1;
            else
                t = (1+sqrt(1+4*obj.tOld^2))/2;
                beta = obj.tOldOld/t;  % beta = i/(i+3);
                obj.tOldOld = obj.tOld;
                obj.tOld = t;
            end
        end
        
        function computeGradientStep(obj)
            tauV = obj.tau;
            D = obj.Differential;
            g = obj.uStar;
            p = obj.pDual;
            p = p - tauV*D*(D'*p - g);
            obj.pDual = p;
        end
        
        
        function projectInTheBall(obj)
            lam = obj.lambda;
            normP = obj.computeNormP();
            no  = max(1,normP/lam);
            p = obj.pDual;
            p = p./[no;no];
            obj.pDual = p;
        end
        
        function normP = computeNormP(obj)
            p = obj.pDual;
            mn = obj.rowsColumns;
            normP = hypot(p(1:mn),p(mn+1:end));
        end
        
        function computeU(obj)
            g = obj.uStar;
            D = obj.Differential;
            p = obj.pDual;
            obj.u = g(:) - D'*p;
        end
        

        
        function printIterationImage(obj,i)
            if mod(i,5)==0
                display(i)
                subFigNum = 2;
                obj.plotImage(obj.u,subFigNum)
            end
        end
        
        function E = computeEnergy(obj)
            p = obj.pDual;
            D = obj.Differential;
            g = obj.uStar;
            gN = obj.gNorm;
            r = D'*p - g;
            E = 0.5*(r'*r)/gN;
        end
        
        function gap = computeDualGap(obj)
            lam = obj.lambda;
            ut = obj.u;
            p = obj.pDual;
            D = obj.Differential;
            gN = obj.gNorm;            
            g = obj.uStar;
            q = D'*p;
            gap = lam*sum(abs(D*ut)) + q'*ut;
            gap = gap/gN;
        end
        
        function plotInitialImage(obj)
            obj.plotImage(obj.u0,1);
            obj.plotImage(obj.uStar,3);
        end
        
        function plotImage(obj,u,subFigNum)
            figure(1)
            m = obj.rows;
            n = obj.columns;
            figure(1);
            subplot(1,3,subFigNum)
            up = reshape(u,m,n);
            imagesc(up);
            colormap(gray); drawnow();
        end
        
        function plotEnergy(obj)
            figure(2)
            plot(1:obj.maxIter,[obj.En(:),obj.EnA(:)],'-+')
            legend('Forward Backward','Accelerated ForwardBackward')
        end
        
        function plotDualGap(obj)
            figure(3)
            plot(1:obj.maxIter,[obj.gap(:),obj.gapA(:)],'-+')
            legend('Forward Backward','Accelerated ForwardBackward')
        end

    end

    
end