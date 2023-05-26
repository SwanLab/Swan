classdef CorrectorCoefficientsComputer < handle
    
    properties (Access = private)
        quadrature
        dVolum
        oCorrectorDerivative        
    end
    
    properties (Access = private)
       mesh
       orthogonalCorrector       
    end
    
    methods (Access = public)
        
        function obj = CorrectorCoefficientsComputer(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createDvolum();
        end

        function c = compute(obj,b)
            obj.createOrthogonalCorrectorDerivatives();          
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(b);            
           c = LHS\RHS;
        end                    
                   
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh                = cParams.mesh;
            obj.orthogonalCorrector = cParams.orthogonalCorrector;
        end
        
        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            obj.quadrature = q;
        end

        function createOrthogonalCorrectorDerivatives(obj)
            psi = obj.orthogonalCorrector;
            q   = obj.quadrature;
            nSing = numel(psi);
            dP = cell(nSing,1);
            for iSing = 1:nSing
                dPsiV = psi{iSing}.computeGradient(q);
                dP{iSing} = dPsiV.fValues;
            end
            obj.oCorrectorDerivative = dP;
        end
        
        function createDvolum(obj)
            q = obj.quadrature;
            nDim = obj.mesh.ndim;
            dV(1,:,:) = obj.mesh.computeDvolume(q);
            dV = repmat(dV,nDim,1,1);
            obj.dVolum = dV;
        end

        function LHS = computeLHS(obj)
            dP  = obj.oCorrectorDerivative;
            dV  = obj.dVolum;
            nSing = numel(obj.orthogonalCorrector);            
            LHS = zeros(nSing,nSing);            
            for iS = 1:nSing
                dPi = dP{iS};
                for jS = 1:nSing
                    dPj = dP{jS};                    
                    lhs   = dPi.*dPj.*dV;
                    LHS(iS,jS) = sum(lhs(:));
                end
            end
        end

        function RHS = computeRHS(obj,b)
           bG    = obj.computeOrientationInGauss(b);
           nSing = numel(obj.orthogonalCorrector);              
           RHS = zeros(nSing,1);    
           dP = obj.oCorrectorDerivative;
           dV   = obj.dVolum;           
           for iS = 1:nSing
                dPi = dP{iS};               
                rhs = dPi.*bG.*dV;
                RHS(iS) = sum(rhs(:));
            end  
        end

        function bfG = computeOrientationInGauss(obj,b)
            q      = obj.quadrature;
            xGauss = q.posgp;
            bfG    = b.evaluate(xGauss);            
            %bfG    = permute(b.evaluate(xGauss),[1 3 2]);
        end
 
    end
    
end