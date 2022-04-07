classdef Constraint < handle
    
    properties (GetAccess = public, SetAccess = private)
        gradient
        value
    end

    properties (Access = private) % classes
        designVariable
        bendingMat
        stiffnessMat
        eigModes
    end
    
    properties (Access = private) % computed
        bendingMatrix
        elementalBendingMatrix
        stiffnessMatrix         
        v1
        v2
        V
        lambda
        D
    end
    
    properties (Access = private) % input  
       nElem
       freeNodes
       length
       youngModulus
       inertiaMoment
       nConstraints
    end
    
    methods (Access = public)

        function obj = Constraint(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj,iter)
             obj.computeStiffnessMatrix();
             obj.computeBendingMatrix();
             obj.computeEigModes(iter);           
            obj.computeFunctions();
            obj.computeGradients();
        end

        function plotModes(obj)
            x = obj.designVariable.value;
            obj.eigModes.plot(x);
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.nElem          = cParams.nElem;
            obj.freeNodes      = cParams.freeNodes;
            obj.length         = cParams.length;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;   
            obj.nConstraints   = cParams.nConstraints;
            obj.bendingMat     = cParams.settings.bendingMat;
            obj.stiffnessMat   = cParams.settings.stiffnessMat;
            obj.eigModes       = cParams.settings.eigMod;
        end

        function computeStiffnessMatrix(obj)
            obj.stiffnessMat.compute();
            obj.stiffnessMatrix =  obj.stiffnessMat.stiffnessMatrix;
        end

        function computeBendingMatrix(obj)
            obj.bendingMat.compute();
            obj.bendingMatrix =  obj.bendingMat.bendingMatrix;
            obj.elementalBendingMatrix = obj.bendingMat.elementalBendingMatrix;
        end

        function computeEigModes(obj,iter)
            K = obj.stiffnessMatrix;
            B = obj.bendingMatrix;
            obj.eigModes.compute(K,B,iter);
            obj.lambda = obj.eigModes.lambda; 
            obj.D      = obj.eigModes.D;
            obj.V      = obj.eigModes.V;
            obj.v1     = obj.eigModes.v1;
            obj.v2     = obj.eigModes.v2;
        end

        function fx = computeFunctions(obj) 
            x = obj.designVariable.value;
            l = obj.lambda;
            N = obj.nElem;
            fx = [x(N+1)-l(1),x(N+1)-l(2),(1/N)*sum(x(1:N))-1]';
            obj.value = fx;
        end

        function [dfdx, dfdx2] = computeGradients(obj)
                Belem = obj.elementalBendingMatrix;
                x = obj.designVariable.value;
                N = obj.nElem;
                m = obj.nConstraints;
                dfdx=zeros(m,N+1);
                dfdx(3,1:N)=(1/N)*ones(1,N);
                dfdx2 = 0*dfdx;
                if abs(obj.D(2,2)-obj.D(1,1))> 1 % Simple eigenvalues
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
                else % Dobles eigenvalues
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
                    % DERIVATIVES MATRIX DEFINITION
                    for i=1:N
                        dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        A=[dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                        [U,R]=eigs(A,2,'SM');
                        S=sort(diag(R));
                        dfdx(1,i)=-S(1);
                        dfdx(2,i)=-S(2);
                    end
                end  
                dfdx(1,N+1)=1;
                dfdx(2,N+1)=1; 
                obj.gradient = dfdx';
        end
    end
    
end