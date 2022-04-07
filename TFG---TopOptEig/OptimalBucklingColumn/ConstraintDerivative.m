classdef ConstraintDerivative < handle

    properties (GetAccess = public, SetAccess = private)
        dfdx
        dfdx2
    end
    
    properties (Access = private)
       nElem
       nConstraint
    end
    
    methods (Access = public)
        
        function obj = ConstraintDerivative(cParams)
            obj.init(cParams)
        end
        
        function computeSimple(obj,x,Be,v1,v2)
            obj.simple(x,Be,v1,v2);
        end
        
        function computeDouble(obj,x,Be,V)
            obj.double(x,Be,V);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem       = cParams.nElem;
            obj.nConstraint = cParams.m;
            obj.dfdx = zeros(obj.nConstraint,obj.nElem+1);
            obj.dfdx(3,1:obj.nElem)=(1/obj.nElem)*ones(1,obj.nElem);
            obj.dfdx2 = 0*obj.dfdx;
        end
        
        function simple(obj,x,Belem,v1,v2)
            N = obj.nElem;
            W=zeros(2*N+2,2);
            for i=3:2*N
                W(i,1)=v1(i-2);
            end
            for i=1:N
                obj.dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,1));
            end
            for i=3:2*N
                W(i,2)=v2(i-2);
            end
            for i=1:N
                obj.dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
            end
        end
        
        function double(obj,x,Belem,V)
            N = obj.nElem;
            Q1=zeros(2*N+2,1);
            Q2=zeros(2*N+2,1);
            dQ1=zeros(N,1);
            dQ2=zeros(N,1);
            dQ1Q2=zeros(N,1);
            for i=3:2*N
                Q1(i,1) = V(i-2,1);
            end
            for i=3:2*N
                Q2(i,1) = V(i-2,2);
            end
            % DERIVATIVES MATRIX DEFINITION
            for i=1:N
                dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                A=[dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                [U,R]=eigs(A,2,'SM');
                S=sort(diag(R));
                obj.dfdx(1,i)=-S(1);
                obj.dfdx(2,i)=-S(2);
            end
        end
       
        
    end
    
end