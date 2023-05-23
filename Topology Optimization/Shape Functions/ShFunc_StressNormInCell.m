classdef ShFunc_StressNormInCell < ShFunWithElasticPdes
    
    % Does not work anymore
    properties (Access = private)
        pNorm = 2;
        stressHomog
        strainHomog
        Chomog
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNormInCell(settings)
            obj@ShFunWithElasticPdes(settings);
            obj.createEquilibriumProblem(settings.filename);
            obj.stressHomog = settings.stressHomog;
        end
        
        function v = getValue(obj)
            v = obj.value;
        end
        
        function setVstrain(obj,s)
            sV(1,:) = s;
            obj.physProb.element.setVstrain(sV);
        end
        
        function setPnorm(obj,p)
            obj.pNorm = p;
        end
        function computeCostAndGradient(obj,x)
            obj.updateMaterialProperties(x);
            obj.solvePDEs();
            obj.computeFunctionValue();
            obj.computeGradient();
        end
        
        function c = computeCostWithFullDomain(obj)
            nnodes = size(obj.physProb.mesh.coord,1);
            x = -ones(nnodes,1);
            obj.computeCostAndGradient(x);
            c = obj.value;
        end
        
        function s = computeMaxStressWithFullDomain(obj)
            nnodes = size(obj.physProb.mesh.coord,1);
            x = -ones(nnodes,1);
            obj.updateMaterialProperties(x);
            obj.solveCellProblem();
            stress = obj.physProb.variables.stress;
            ngaus = obj.physProb.element.quadrature.ngaus;
            nstre = obj.physProb.element.getNstre();
            V = sum(sum(obj.physProb.geometry.dvolu));
            dV = obj.physProb.element.geometry.dvolu;
            s = obj.obtainMaxSigmaNorm(stress,ngaus,nstre,dV,V);
        end
        
    end
    
    methods (Access = protected)
        
        function updateGradient(obj)
        end
        
        function computeGradient(obj)
            obj.gradient = 0;
        end
        
        function solvePDEs(obj)
           obj.computeHomogenizedTensor();
           obj.computeHomogenizedStrain();
           obj.solveCellProblem();
        end
        
        function computeFunctionValue(obj)
            stress       = obj.physProb.variables.stress;
            stress_fluct = obj.physProb.variables.stress_fluct;
            v = obj.integrateStressNorm(obj.physProb);
            obj.value = v;
        end
        
    end
    
    methods (Access = private)
        
        function solveCellProblem(obj)
            cellProblem = obj.physProb;
            sH(1,:) = obj.strainHomog;
            cellProblem.element.setVstrain(sH);
            cellProblem.setMatProps(obj.matProps);
            cellProblem.computeVariables();
        end
        
        function computeHomogenizedTensor(obj)
            cellProblem = obj.physProb;
            Ch = cellProblem.solve();
            obj.Chomog = Ch;
        end
        
        function computeHomogenizedStrain(obj)
             Ch = obj.Chomog;
             stress = obj.stressHomog;
             strain = Ch\stress;
             obj.strainHomog = strain;
        end
        
        function value = integrateStressNorm(obj,physProb)
            nstre = physProb.element.getNstre();
            V = sum(sum(physProb.geometry.dvolu));
            ngaus = physProb.element.quadrature.ngaus;
            dV    = physProb.element.geometry.dvolu;
            stress = obj.physProb.variables.stress;
            value = obj.integratePNormOfL2Norm(stress,ngaus,nstre,dV,V);
        end
        
        
        function v = integratePNormOfL2Norm(obj,stress,ngaus,nstre,dV,V)
            p = obj.pNorm;
            v = 0;
            for igaus = 1:ngaus
                s = zeros(size(stress,3),1);
                for istre = 1:nstre
                    Si = squeeze(stress(igaus,istre,:));
                    factor = obj.computeVoigtFactor(istre,nstre);
                    Sistre = factor*(Si.^2);
                    s = s + Sistre;
                end
                s = sqrt(s);
                sigmaNorm = s.^p;
                v = v + 1/V*sigmaNorm'*dV(:,igaus);
            end
        end
        
        function v = obtainMaxSigmaNorm(obj,stress,ngaus,nstre,dV,V)
            v = [];
            for igaus = 1:ngaus
                s = zeros(size(stress,3),1);
                for istre = 1:nstre
                    Si = squeeze(stress(igaus,istre,:));
                    factor = obj.computeVoigtFactor(istre,nstre);
                    Sistre = factor*(Si.^2);
                    s = s + Sistre;
                end
                v = max([sqrt(s);v]);
            end
            
        end
        
        function v = integratePNormOfComponents(obj,stress,ngaus,nstre,dV,V)
            p = obj.pNorm;
            v = 0;
            for igaus = 1:ngaus
                s = zeros(size(stress,1),size(stress,3));
                for istre = 1:nstre
                    Si = squeeze(stress(igaus,istre,:));
                    factor = obj.computeVoigtFactor(istre,nstre);
                    Sistre = factor*(Si.^p);
                    s = s + Sistre;
                end
                v = v + 1/V*s'*dV(:,igaus);
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