classdef ShFunc_InternalEnergy < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        materialPhaseField
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergy(cParams)
            obj.init(cParams)            
        end
        
        function E = computeFunction(obj,u,phi)

            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            
            e = u.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            
            s.quadrature = quad;
            s.phi = phi;
            s.derivative = deriv;
            obj.materialPhaseField.computeMatInt(s);

            q.mesh = obj.mesh;
            q.type = 'InternalEnergy';
            q.quadType = 'QUADRATIC';
            int = Integrator.create(q);

            C = obj.materialPhaseField.material.C;
            e = u.computeSymmetricGradient(quad);
            E = 0.25*int.compute(e,C);
        end
        
        function J = computeGradient(obj,u,phi)
            DenergyFun =  obj.createAbstractDerivativeEnergyFunction('LINEAR',1);
            test = P1Function.create(obj.mesh,1);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 'LINEAR';
            RHS = RHSintegrator.create(s);
            obj.Fi = RHS.compute(DenergyFun,test);
            
        end
        
        function H = computeHessian(obj,u,phi)
            DDenergyFun =  obj.createAbstractDerivativeEnergyFunction('LINEAR',2);
            
            s.function = DDenergyFun;
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            LHS = LHSintegrator.create(s);
            obj.Mi = LHS.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function energyFun = createAbstractDerivativeEnergyFunction(obj,quadOrder,deriv,u,phi)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            
            e = u.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            
            s.quadrature = quad;
            s.phi = phi;
            s.derivative = deriv;
            obj.materialPhaseField.computeMatInt(s);
            C = obj.materialPhaseField.material.C;
            DDenergyVal = obj.computeEnergyField(e,C);
            
            s.fValues = DDenergyVal;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            energyFun = FGaussDiscontinuousFunction(s);
        end
        
        function energyVal = computeEnergyField(obj,e,C)
            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            energyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eI = squeeze(e.fValues(iStre,iGauss,:));
                        eJ = squeeze(e.fValues(jStre,iGauss,:));
                        Cij = squeeze(C(iStre,jStre,:,iGauss));
                        eStre(1,1,:) = (eI.*Cij.*eJ)';
                        energyVal(1,iGauss,:) = energyVal(1,iGauss,:) + eStre;
                    end
                end
            end
            energyVal = 0.5*energyVal;
        end
        
        
    end
    
end