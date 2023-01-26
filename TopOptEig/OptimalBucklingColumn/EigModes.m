classdef EigModes < handle
       
    properties (Access = private)
        mesh
        V
        v1
        v2
        D
        eigModesPlotter
        lambda
        mode1Disp
        mode2Disp
        freeNodes        
    end

    properties (Access = private)
        dim
        sectionVariables
        stiffnessMatComputer
        bendingMatComputer
        inertiaMoment
        youngModulus
    end

    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
            obj.createDimensions();
            obj.createBoundaryConditions();
            obj.createStiffnessMatrix();
            obj.createBendingMatrix();
            obj.createEigModesPlotter();
        end             

        function plot(obj,A,iter)
            p = obj.eigModesPlotter;
            [m1,m2] = obj.computeBucklingModes(obj.v1,obj.v2);
            obj.mode1Disp = m1;
            obj.mode2Disp = m2;
            p.plot(A,m1,m2,iter,obj.D)
        end

        function l = provideEigenValue(obj)
            obj.computeEigenModesAndValues();            
            obj.computeLambda();   
            l = obj.lambda;
%             gamma = obj.designVariable.getFirstEigenMode();            
%             fx = gamma-obj.lambda(eigNum);
        end

       function grad = provideDerivative(obj,eigNum)
            obj.reorderModes(obj.lambda,obj.V,obj.D);
            Belem =  obj.bendingMatComputer.elementalBendingMatrix;
            %x = obj.designVariable.getColumnArea();
            %x = obj.designVariable.getColumnRadius();
            nElem = obj.mesh.nelem;
            eigV1 = obj.D(1,1);
            eigV2 = obj.D(2,2);
            difEigs = abs(eigV2-eigV1);
            if difEigs > 1 
                dfdx = obj.computeSimpleEig(Belem);
            else 
                dfdx = obj.computeDoubleEig(Belem);
            end
            dfdx(1,nElem+1) = 1;
            dfdx(2,nElem+1) = 1;
            grad = dfdx(eigNum,:);
        end

    end

    methods (Access = private)
        
        function dfdx = computeSimpleEig(obj,Belem)
            d = obj.dim;
            free = obj.freeNodes;
            ndofe = d.ndofsElem;
            ndofn = d.ndofsElem/d.nnodeElem;
            W = zeros(d.ndofs,2);
            W(free,1) = obj.v1;
            W(free,2) = obj.v2;
            nElem = obj.mesh.nelem;
            dI = obj.sectionVariables.computeInertiaDerivative();
            for i = 1:nElem
                index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
                dx = dI(i,1);
                %dx = 2*x(i,1);
                %dx = 4*pi^2*x(i,1).^3;
                dfdx(1,i) = -dx*(W(index,1)'*Belem(:,:,i)*W(index,1));
                dfdx(2,i) = -dx*(W(index,2)'*Belem(:,:,i)*W(index,2));
            end
%             for i = 1:nElem
%                 index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
%                 Eig1(:,nElem) = W(index,1);
%                 Eig2(:,nElem) = W(index,2);
%             end
%             for i = 1:ndofe
%                 for j = 1:ndofe
%                     Bij = squeeze(Belem(i,j,:));
%                     Eig1j = Eig1(j,:)';
% 
%                 end
%             end

        end

        function dfdx = computeDoubleEig(obj,Belem)
            d    = obj.dim;
            free = obj.freeNodes;
            ndofe = d.ndofsElem;
            ndofn = d.ndofsElem/d.nnodeElem;
            nElem = obj.mesh.nelem;
            W1    = zeros(d.ndofs,1);
            W2    = zeros(d.ndofs,1);
            dW1   = zeros(nElem,1);
            dW2   = zeros(nElem,1);
            dW1W2 = zeros(nElem,1);
            W1(free,1) = obj.v1;
            W2(free,1) = obj.v2;
            dI = obj.sectionVariables.computeInertiaDerivative();
            x = obj.sectionVariables.computeArea();
            for i=1:nElem
                index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
                dx = dI(i,1);
                %dx = 2*x(i,1);
                %dx = 4*pi^2*x(i,1).^3;
                dW1(i,1)= dx*(W1(index,1)'*Belem(:,:,i)*W1(index,1));
                dW2(i,1)= dx*(W2(index,1)'*Belem(:,:,i)*W2(index,1));
                dW1W2(i,1)= (2*x(i,1))*(W1(index,1)'*Belem(:,:,i)*W2(index,1));
                A = [dW1(i,1) dW1W2(i,1); dW1W2(i,1) dW2(i,1)];
                A1 = dW1(i,1);
                A12 = dW1W2(i,1);
                A21 = dW1W2(i,1);
                A2 = dW2(i,1);
                a=1;
                b = -(A1 + A2);
                c = A1*A2 - A12*A21;
                lambd1 = -b+sqrt(b^2-4*a*c)/(2*a);
                lambd2 = -b-sqrt(b^2-4*a*c)/(2*a);
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
            %obj.designVariable       = cParams.designVariable;
            obj.sectionVariables    = cParams.sectionVariables;
            obj.inertiaMoment        = cParams.inertiaMoment;
            obj.youngModulus         = cParams.youngModulus;
        end
        
        function createDimensions(obj)
            s.mesh = obj.mesh;
            s.pdim = '2D';
            s.ngaus = 2;
            s.type = 'Vector';
            s.name = 'x';
            s.ndimf = 2;
            s.fieldName = 'u';
            d = DimensionVariables.create(s);
            d.compute();
            obj.dim = d;
        end
        
         function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.ndofs-1,d.ndofs]);
            %fixnodes = [1,d.ndofs-1];
            nodes = 1:d.ndofs;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
         end

        function createStiffnessMatrix(obj)
            s.type = 'StiffnessMatrixColumn';
            s.dim = obj.dim;
            s.mesh = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.freeNodes      = obj.freeNodes;
            K = LHSintegrator.create(s);
            obj.stiffnessMatComputer = K;
        end

        function createBendingMatrix(obj)
            s.type         = 'BendingMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            %s.designVariable = obj.designVariable;
            s.sectionVariables = obj.sectionVariables;
            s.freeNodes      = obj.freeNodes;
            B = LHSintegrator.create(s);
            obj.bendingMatComputer = B;
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

        function computeLambda(obj)
            l = sort(diag(obj.D));
            obj.lambda = l;
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