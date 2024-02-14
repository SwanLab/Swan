classdef ShFunc_InternalEnergy < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        materialPhaseField
        stressFun
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergy(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quad)
            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad);
            energyFun =  obj.createEnergyFunction(quad,u);
            s.mesh = obj.mesh;
            s.type = 'Function';
            s.quadType = quad.order;
            int = Integrator.create(s);
            F = 0.5*int.compute(energyFun);
        end
        
        function [Ju, Jphi] = computeGradient(obj,u,phi,quad)
            obj.materialPhaseField.computeFirstDerivativeInterpolatedMaterial(phi,quad); 
            dEnergyFun =  obj.createEnergyFunction(quad,u);
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quad.order;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);

            
            sigma = obj.createStressFunction(u,phi,quad);
            test = LagrangianFunction.create(obj.mesh, 2, 'P1');

            s.mesh = obj.mesh;
            s.quadratureOrder = 'LINEAR';
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end
        
        function [Huu, Hphiphi] = computeHessian(obj,u,phi,quad)
            obj.materialPhaseField.computeSecondDerivativeInterpolatedMaterial(phi,quad); 
            ddEnergyFun =  obj.createEnergyFunction(quad,u);
            
            s.function = ddEnergyFun;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            Hphiphi = 0.5*LHS.compute();

            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = obj.materialPhaseField.material;
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
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

        function stressFun = createStressFunction(obj,u,phi,quad)
            s.mesh = obj.mesh;
            s.materialPhaseField = obj.materialPhaseField;
            s.u = u;
            s.phi = phi;
            s.quadrature = quad;
            stressFun = StressFunctions(s);
        end
    end
    
end