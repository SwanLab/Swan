classdef EigModes < handle
    
    properties (Access = public)
         
    end
    
    properties (Access = private)
        length


        V
        v1
        v2
        D
        eigModesPlotter
        lambda
    end

    properties (Access = private)
        designVariable
        stiffnessMatComputer
        bendingMatComputer
    end

    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
            obj.createEigModesPlotter();
        end             

        function plot(obj,A,iter)
            p = obj.eigModesPlotter;
            [m1,m2] = obj.computeBucklingModes(obj.v1,obj.v2);
            p.plot(A,m1,m2,iter,obj.D)
        end

        function fx = provideFunction(obj,eigNum)
            obj.computeEigenModesAndValues();            
            obj.lambda = obj.computeLambda();                
            gamma = obj.designVariable.getFirstEigenMode();            
            fx = gamma-obj.lambda(eigNum);
        end

        function grad = provideDerivative(obj,eigNum)
            obj.reorderModes(obj.lambda,obj.V,obj.D);
            Belem =  obj.bendingMatComputer.elementalBendingMatrix;
            x = obj.designVariable.getColumnArea;
            N = obj.designVariable.getNelem();
            if abs(obj.D(2,2)-obj.D(1,1))> 1
                W=zeros(2*N+2,2);
                for i=3:2*N
                    W(i,1)=obj.v1(i-2);
                end
                for i=1:N
                    dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,1));
                end
                for i=3:2*N
                    W(i,2)=obj.v2(i-2);
                end
                for i=1:N
                    dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
                end
            else
                obj.D
                disp('dobles')
                Q1=zeros(2*N+2,1);
                Q2=zeros(2*N+2,1);
                dQ1=zeros(N,1);
                dQ2=zeros(N,1);
                dQ1Q2=zeros(N,1);
                for i=3:2*N
                    Q1(i,1)=obj.V(i-2,1);
                end
                for i=3:2*N
                    Q2(i,1)=obj.V(i-2,2);
                end
                A  = zeros(2,2);
                for i=1:N
                    dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                    dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                    dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                    A = [dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                    [U,R] = eigs(A,2,'SM');
                    S = sort(diag(R));
                    dfdx(1,i)=-S(1);
                    dfdx(2,i)=-S(2);
                end
            end
            dfdx(1,N+1)=1;
            dfdx(2,N+1)=1;
            grad = dfdx(eigNum,:);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.length         = cParams.length;
            obj.designVariable = cParams.designVariable;
            obj.stiffnessMatComputer = cParams.stiffnessMatComputer;
            obj.bendingMatComputer   = cParams.bendingMatComputer;
        end

        function computeEigenModesAndValues(obj) 
            obj.stiffnessMatComputer.compute();
            obj.bendingMatComputer.compute();            
            Kfree  = obj.stiffnessMatComputer.provideFreeStiffnessMatrix();
            Bfree  = obj.bendingMatComputer.provideFreeBendingMatrix();
            obj.computeEigenFunctionAndValues(Bfree,Kfree);         
        end

        function createEigModesPlotter(obj)
            s.length = obj.length;
            p = EigModesPlotter(s);
            obj.eigModesPlotter = p;
        end

        function l = computeLambda(obj)
            l = sort(diag(obj.D));       
        end

        function computeEigenFunctionAndValues(obj,B,K)
            [v,d] = eigs(B,K,2,'SM');
            obj.V  = v;
            obj.D  = d; 
        end
        
        function reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                v1=V(:,1);
                v2=V(:,2);
            else
                v1=V(:,2);
                v2=V(:,1);
            end
            obj.v1 = v1;
            obj.v2 = v2;             
        end

        function [m1,m2] = computeBucklingModes(obj,v1,v2)
            N = obj.designVariable.getNelem();
            Mode1=zeros(2*(N+1),1);
            Mode2=zeros(2*(N+1),1);
            for i=3:2*N
                Mode1(i)=v1(i-2);
                Mode2(i)=v2(i-2);
            end
            m1 = Mode1;
            m2 = Mode2;
        end             


          

    end
    
end