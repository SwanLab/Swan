classdef ShFunc_StressNorm < ShFunWithElasticPdes
    
    properties (Access = private)
        pNorm = 2;
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm(settings)
            obj@ShFunWithElasticPdes(settings);
            obj.createEquilibriumProblem(settings.filename);
        end
        
        function v = getValue(obj)
            v = obj.value;
        end
        
        function setVstrain(obj,s)
            sV(1,:) = s;
            obj.physProb.element.setVstrain(sV);
        end
        function computeCostAndGradient(obj,x)
            obj.updateMaterialProperties(x);
            obj.solvePDEs();
            obj.computeFunctionValue();
            obj.computeGradient();
        end
        
    end
    
    methods (Access = protected)
        
        function updateGradient(obj)
        end
        
        function computeGradient(obj)
            obj.gradient = 0;
        end
        
        function solvePDEs(obj)
            obj.physProb.setMatProps(obj.matProps);
            obj.physProb.computeVariables();
        end
        
        function computeFunctionValue(obj)
            stress       = obj.physProb.variables.stress;
            stress_fluct = obj.physProb.variables.stress_fluct;
            v = obj.integrateStress(obj.physProb);
            obj.value = v;
        end
        
    end
    
    methods (Access = private)
        
        function v = integrateStress(obj,physProb)
            nstre = physProb.element.getNstre();
            V = sum(sum(physProb.geometry.dvolu));
            ngaus = physProb.element.quadrature.ngaus;
            dV    = physProb.element.geometry.dvolu;
            
            stress = obj.physProb.variables.stress;
            
            p = obj.pNorm;
            v = 0;
            for istre = 1:nstre
                for igaus = 1:ngaus
                    Si = squeeze(stress(igaus,istre,:));
                    factor = obj.computeVoigtFactor(istre,nstre);
                    vij = 1/V*factor*(Si.^p)'*dV(:,igaus);
                    v = v + vij;
                end
            end
            
        end
        
        function f = computeVoigtFactor(obj,k,nstre)
            if nstre == 3
                if k < 3
                    f = 1;
                else
                    f = 2;
                end
            elseif nstre ==6
                if k <= 3
                    f = 1;
                else
                    f = 2;
                end
            end
        end
        
        
    end
    
end