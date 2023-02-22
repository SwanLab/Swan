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
            eigV1 = obj.D(1,1);
            eigV2 = obj.D(2,2);
            difEigs = abs(eigV2-eigV1);
            if difEigs > 1 
                dfdx = obj.computeSimpleEig(Belem);
            else 
                dfdx = obj.computeDoubleEig(Belem);
            end
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
            
            nVar = obj.sectionVariables.nDesVarElem;
            gElemt = ndofn*((1:nElem)-1);
            dx=dI;
            WB = zeros(2,ndofe,nElem);
            WBW = zeros(nElem,2);
            for iEig = 1:2
                for iDof = 1:ndofe
                    for jDof = 1:ndofe
                        indexJ = gElemt(:)+jDof;
                        Bij(:,1) = squeeze(Belem(jDof,iDof,:));
                        WBj(1,1,:) = W(indexJ,iEig).*Bij;
                        WB(iEig,iDof,:) = WB(iEig,iDof,:)+WBj;
                    end
                    indexI = gElemt(:)+iDof;
                    auxA = squeeze(WB(iEig,iDof,:));
                    WBW(:,iEig) = auxA.*W(indexI,iEig)+WBW(:,iEig);
                end
                WBWt = repmat(WBW(:,iEig),nVar,1);
                dfdx(iEig,:) = -dx.*WBWt;
            end
            dfdx(1,nVar*nElem+1) = 1;
            dfdx(2,nVar*nElem+1) = 1;
        end


        function dfdx = computeDoubleEig(obj,Belem)
            d    = obj.dim;
            free = obj.freeNodes;
            ndofe = d.ndofsElem;
            ndofn = d.ndofsElem/d.nnodeElem;
            nElem = obj.mesh.nelem;
            W1    = zeros(d.ndofs,1);
            W2    = zeros(d.ndofs,1);
            dW1   = zeros(2*nElem,1);
            dW2   = zeros(2*nElem,1);
            dW1W2 = zeros(2*nElem,1);
            W1(free,1) = obj.v1;
            W2(free,1) = obj.v2;
            dI = obj.sectionVariables.computeInertiaDerivative();
            switch obj.sectionVariables.nDesVarElem
                case 1
                    for iElem=1:nElem
                        index = ndofn*(iElem-1)+1: ndofn*(iElem-1)+ndofe;
                        dx = dI(iElem,1);
                        dW1(iElem,1)= dx*(W1(index,1)'*Belem(:,:,iElem)*W1(index,1));
                        dW2(iElem,1)= dx*(W2(index,1)'*Belem(:,:,iElem)*W2(index,1));
                        dW1W2(iElem,1)= dx*(W1(index,1)'*Belem(:,:,iElem)*W2(index,1));
        
                        A = [dW1(iElem,1) dW1W2(iElem,1); dW1W2(iElem,1) dW2(iElem,1)];
                        A1 = dW1(iElem,1);
                        A12 = dW1W2(iElem,1);
                        A21 = dW1W2(iElem,1);
                        A2 = dW2(iElem,1);
                        a=1;
                        b = -(A1 + A2);
                        c = (A1*A2)-(A12*A21);
                        lambd1 = (-b+sqrt(b^2-4*a*c))/(2*a);
                        lambd2 = (-b-sqrt(b^2-4*a*c))/(2*a);
                        %S = sort([lambd1,lambd2]);
                        [U,R] = eigs(A,2,'SM');
                        S = sort(diag(R));
                        dfdx(1,iElem) = -S(1);
                        dfdx(2,iElem) = -S(2);
                    end
                    dfdx(1,nElem+1) = 1;
                    dfdx(2,nElem+1) = 1;
                case 2
                    for iElem=1:nElem
                        index = ndofn*(iElem-1)+1: ndofn*(iElem-1)+ndofe;
                        da = dI(iElem,1);
                        db = dI(nElem+iElem,1);
                        dW1(iElem,1)= da*(W1(index,1)'*Belem(:,:,iElem)*W1(index,1));
                        dW2(iElem,1)= da*(W2(index,1)'*Belem(:,:,iElem)*W2(index,1));
                        dW1W2(iElem,1)= da*(W1(index,1)'*Belem(:,:,iElem)*W2(index,1));
                        dW1(nElem+iElem,1)= db*(W1(index,1)'*Belem(:,:,iElem)*W1(index,1));
                        dW2(nElem+iElem,1)= db*(W2(index,1)'*Belem(:,:,iElem)*W2(index,1));
                        dW1W2(nElem+iElem,1)= db*(W1(index,1)'*Belem(:,:,iElem)*W2(index,1));
                    end
                    
                    for iElem=1:2*nElem
                        A = [dW1(iElem,1) dW1W2(iElem,1); dW1W2(iElem,1) dW2(iElem,1)];
                        A1 = dW1(iElem,1);
                        A12 = dW1W2(iElem,1);
                        A21 = dW1W2(iElem,1);
                        A2 = dW2(iElem,1);
                        a=1;
                        b = -(A1 + A2);
                        c = A1*A2 - A12*A21;
                        lambd1 = -b+sqrt(b^2-4*a*c)/(2*a);
                        lambd2 = -b-sqrt(b^2-4*a*c)/(2*a);
                        [U,R] = eigs(A,2,'SM');
                        S = sort(diag(R));
                        dfdx(1,iElem) = -S(1);
                        dfdx(2,iElem) = -S(2);
                    end
                    dfdx(1,2*nElem+1) = 1;
                    dfdx(2,2*nElem+1) = 1;
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
            obj.stiffnessMatComputer.compute();            
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