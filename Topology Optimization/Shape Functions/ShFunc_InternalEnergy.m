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
        
        function F = computeFunction(obj,u,phi,quad)
            e = u.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad);
            C = obj.materialPhaseField.material.C;

            s.mesh = obj.mesh;
            s.type = 'InternalEnergy';
            s.quadType = quad.order;
            int = Integrator.create(s);
            F = 0.5*int.compute(e,C);

            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad); %% Redundant
            energyFun =  obj.createEnergyFunction(quad,u);
            s.mesh = obj.mesh;
            s.type = 'Function';
            s.quadType = quad.order;
            int2 = Integrator.create(s);
            F2 = 0.5*int2.compute(energyFun);
        end
        
        function J = computeGradient(obj,u,phi,quad)
            obj.materialPhaseField.computeFirstDerivativeInterpolatedMaterial(phi,quad); 
            dEnergyFun =  obj.createEnergyFunction(quad,u);
            test = P1Function.create(obj.mesh,1);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quad.order;
            RHS = RHSintegrator.create(s);
            J.phi = 0.5*RHS.compute(dEnergyFun,test);

            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = obj.materialPhaseField.material;
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            J.u = LHS.compute();

        end
        
        function H = computeHessian(obj,u,phi,quad)
            obj.materialPhaseField.computeSecondDerivativeInterpolatedMaterial(phi,quad); 
            ddEnergyFun =  obj.createEnergyFunction(quad,u);
            
            s.function = ddEnergyFun;
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            H.phiphi = 0.5*LHS.compute();

            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = obj.materialPhaseField.material;
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            H.uu = LHS.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.materialPhaseField;
        end
        
        function energyFun = createEnergyFunction(obj,quad,u)         
            e = u.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            
            C = obj.materialPhaseField.material.C;
            energyVal = obj.computeEnergyField(e,C);
            
            s.fValues = energyVal;
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
        end
    end
    
end