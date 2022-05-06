classdef EigModes < handle
    
    properties (Access = public)
         
    end
    
    properties (Access = private)
        mesh
        dim
        V
        v1
        v2
        D
        eigModesPlotter
        lambda
        freeNodes
        mode1Disp
        mode2Disp
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
            obj.mode1Disp = m1;
            obj.mode2Disp = m2;
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
            x = obj.designVariable.getColumnArea();
            nElem = obj.mesh.nelem;
            eigV1 = obj.D(1,1);
            eigV2 = obj.D(2,2);
            difEigs = abs(eigV2-eigV1);
            if difEigs > 1 
                dfdx = obj.computeSimpleEig(Belem,x);
            else 
                dfdx = obj.computeDoubleEig(Belem,x);
            end
            dfdx(1,nElem+1) = 1;
            dfdx(2,nElem+1) = 1;
            grad = dfdx(eigNum,:);
        end

    end

    methods (Access = private)

        function dfdx = computeSimpleEig(obj,Belem,x)
                d = obj.dim;
                free = obj.freeNodes;
                W = zeros(d.ndof,2);
                W(free,1) = obj.v1;
                W(free,2) = obj.v2;
                nElem = obj.mesh.nelem;
                for i = 1:nElem
                    index = 2*(i-1)+1: 2*(i-1)+4;
                    dfdx(1,i) = -(2*x(i,1))*(W(index,1)'*Belem(:,:,i)*W(index,1));
                    dfdx(2,i) = -(2*x(i,1))*(W(index,2)'*Belem(:,:,i)*W(index,2));
                end
        end

        function dfdx = computeDoubleEig(obj,Belem,x)
            d    = obj.dim;
            free = obj.freeNodes;
            ndofe = d.ndofPerElement;
            ndofn = d.ndimField;
            nElem = obj.mesh.nelem;
            Q1    = zeros(d.ndof,1);
            Q2    = zeros(d.ndof,1);
            dQ1   = zeros(nElem,1);
            dQ2   = zeros(nElem,1);
            dQ1Q2 = zeros(nElem,1);
            Q1(free,1) = obj.v1;
            Q2(free,1) = obj.v2;
            for i=1:nElem
                index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
                der_x = 2*x(i,1);
                dQ1(i,1)= der_x*(Q1(index,1)'*Belem(:,:,i)*Q1(index,1));
                dQ2(i,1)= der_x*(Q2(index,1)'*Belem(:,:,i)*Q2(index,1));
                dQ1Q2(i,1)= (2*x(i,1))*(Q1(index,1)'*Belem(:,:,i)*Q2(index,1));
                A = [dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                [U,R] = eigs(A,2,'SM');
                S = sort(diag(R));
                dfdx(1,i) = -S(1);
                dfdx(2,i) = -S(2);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.dim                  = cParams.dim;
            obj.designVariable       = cParams.designVariable;
            obj.stiffnessMatComputer = cParams.stiffnessMatComputer;
            obj.bendingMatComputer   = cParams.bendingMatComputer;
        end

        function computeEigenModesAndValues(obj) 
            obj.stiffnessMatComputer.compute();
            obj.bendingMatComputer.compute();            
            [Kfree,free]  = obj.stiffnessMatComputer.provideFreeStiffnessMatrix();
            obj.freeNodes = free;
            Bfree  = obj.bendingMatComputer.provideFreeBendingMatrix();
            obj.computeEigenFunctionAndValues(Bfree,Kfree);         
        end

        function createEigModesPlotter(obj)
            s.mesh = obj.mesh;
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
                V1=V(:,1);
                V2=V(:,2);
            else
                V1=V(:,2);
                V2=V(:,1);
            end
            obj.v1 = V1;
            obj.v2 = V2;             
        end

        function [m1disp,m2disp] = computeBucklingModes(obj,v1,v2)
            N = obj.mesh.nelem;
            Mode1=zeros(2*(N+1),1);
            Mode2=zeros(2*(N+1),1);
            for i=3:2*N
                Mode1(i)=v1(i-2);
                Mode2(i)=v2(i-2);
            end
            m1 = Mode1;
            m2 = Mode2;
            m1disp = m1(1:2:end);
            m2disp = m2(1:2:end);
        end             

    end
    
end