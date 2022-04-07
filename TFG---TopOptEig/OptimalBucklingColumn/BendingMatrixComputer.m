classdef BendingMatrixComputer < handle
    
    properties (Access = public)
        bendingMatrix
        elementalBendingMatrix 
    end
    
    properties (Access = private)
        nElem
        length
        youngModulus
        inertiaMoment
        designVariable
    end
    
    methods (Access = public)
        
        function obj = BendingMatrixComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeElementalBendingMatrix();
            obj.computeBendingMatrix();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem        = cParams.nElem;
            obj.length       = cParams.length;
            obj.youngModulus = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable;
        end
        
         
        function Be = computeElementalBendingMatrix(obj)
            L = obj.length;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            [c1,c2,c3,c4] = obj.coeffsBending(L,E,I);
            Be(1,1:4) = c1*[c2 c3 -c2 c3];
            Be(2,1:4) = c1*[c3 c4 -c3 c4/2];
            Be(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Be(4,1:4) = c1*[c3 c4/2 -c3 c4];
            obj.elementalBendingMatrix  = Be;
        end
        
        function [c1,c2,c3,c4] = coeffsBending(obj,L,E,I)
            c1 = E*I/L^3;
            c2 = 12;
            c3 = 6*L;
            c4 = 4*L^2;
        end
        
        function computeBendingMatrix(obj)
            x = obj.designVariable.value;
            N = obj.nElem;
            B=sparse(2*N+2, 2*N+2);
            for iElem = 1: N
                iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                B(iDof,iDof)=B(iDof,iDof)+(x(iElem)^2)*obj.elementalBendingMatrix;
            end
            obj.bendingMatrix = B;
        end
    end
    
end