classdef Elastic_Problem_Micro < Elastic_Problem
    %Elastic_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)
    end
    
    %% Private properties definition ======================================
    properties (Access = private)
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem_Micro(problemID)
            obj@Elastic_Problem(problemID);
            obj.dof = DOF_Elastic_Micro(problemID,obj.geometry,obj.mesh);
        end
        
        
        function [Chomog,tstrain,tstress] = computeChomog(obj)
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.interpolation_u.computeShapeDeriv(obj.element.quadrature.posgp)
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
            % obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            nstre = obj.element.getNstre();
            vstrain = diag(ones(nstre,1));
            Chomog =  zeros(nstre,nstre);
            tstrain = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tstress = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            var2print = cell(nstre,1);
            for istre=1:nstre
                obj.element.setVstrain(vstrain(istre,:));
                obj.computeVariables;
                Chomog(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstress(istre,:,:,:) = obj.variables.stress;
                var2print{istre}.stress = obj.variables.stress;
                var2print{istre}.strain = obj.variables.strain;
                var2print{istre}.stress_fluct = obj.variables.strain_fluct;
                var2print{istre}.strain_fluct = obj.variables.strain_fluct;
                var2print{istre}.d_u = obj.variables.d_u;
                var2print{istre}.fext = obj.variables.fext;
            end
            obj.variables.Chomog  = Chomog;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
            
            obj.variables.var2print = var2print;
        end
        
        function P2 = computeAmplificator(obj)
            obj.variables.Ptensor = rand(3,3);
            [Chomog,tstrain,tstres] = obj.computeChomog();
            Ct = StiffnessPlaneStressVoigtTensor;
            Ct.setValue(Chomog);
            St = Inverter.invert(Ct);
            Shomog = St.getValue();
            
            nstre = obj.element.getNstre();
            V = sum(sum(obj.geometry.dvolu));
            ngaus = obj.element.quadrature.ngaus;
            dV    = obj.element.geometry.dvolu;
            C     = obj.element.material.C;
            
            P2  = obj.computeAmplificatorbyPtensors(tstres,V,dV,Shomog,nstre,ngaus);
            P3  = obj.computeAmplificatorsWithStress(tstres,V,dV,Shomog,nstre,ngaus);
            P2b = obj.computeAmplificatorWithAllLoops(tstrain,C,V,dV,Shomog,nstre,ngaus);
        end
        
        
        function P2b = computeAmplificatorWithAllLoops(obj,tstrain,C,V,dV,Shomog,nstre,ngaus)
            P2b = zeros(size(Shomog));
            for istre = 1:nstre
                for jstre = 1:nstre
                    for mstre = 1:nstre
                        for kstre = 1:nstre
                            for lstre = 1:nstre
                                for igaus = 1:ngaus
                                    Pmi = zeros(size(tstrain,4),1);
                                    Pmj = zeros(size(tstrain,4),1);
                                    for qstre = 1:nstre
                                        Cmq = squeeze(C(mstre,qstre,:));
                                        Wkq = squeeze(tstrain(kstre,igaus,qstre,:));
                                        Wlq = squeeze(tstrain(lstre,igaus,qstre,:));
                                        CWmk = Cmq.*Wkq;
                                        CWml = Cmq.*Wlq;
                                        S_ik  = Shomog(kstre,istre);
                                        S_jl  = Shomog(lstre,jstre);
                                        Pmi = Pmi + CWmk*S_ik;
                                        Pmj = Pmj + CWml*S_jl;
                                    end
                                    factor = obj.computeVoigtFactor(mstre,nstre);
                                    P2b(istre,jstre) = P2b(istre,jstre) + 1/V*(Pmi.*Pmj)'*dV(:,igaus)*factor;
                                end
                            end
                            
                        end
                        
                    end
                end
            end
        end
        
        function P = computeAmplificatorbyPtensors(obj,tstres,V,dV,Shomog,nstre,ngaus)
            Pt = obj.computePtensor(tstres,Shomog,nstre);
            P  = obj.computeAmplificatorTensor(Pt,nstre,V,dV,ngaus);
        end
        
        function P = computeAmplificatorTensor(obj,Pt,nstre,V,dV,ngaus)
            P = zeros(nstre,nstre);
            for istre = 1:nstre
                for jstre = 1:nstre
                    for kstre = 1:nstre
                        for igaus = 1:ngaus
                            Pki = squeeze(Pt(kstre,igaus,istre,:));
                            Pkj = squeeze(Pt(kstre,igaus,jstre,:));
                            factor = obj.computeVoigtFactor(kstre,nstre);
                            P(istre,jstre) = P(istre,jstre) + 1/V*factor*(Pki.*Pkj)'*dV(:,igaus);
                        end
                    end
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
        
        
        function Pt = computePtensor(obj,tstres,Shomog,nstre)
            Pt = zeros(size(tstres));
            for istre = 1:nstre
                for jstre = 1:nstre
                    for kstre = 1:nstre
                        Pt(istre,:,jstre,:) = Pt(istre,:,jstre,:) + tstres(istre,:,kstre,:)*Shomog(kstre,jstre);
                    end
                end
            end
        end
        
        
        
        function P = computeAmplificatorsWithStress(obj,tstres,V,dV,Shomog,nstre,ngaus)
            Ch2 = zeros(size(Shomog));
            for istre = 1:nstre
                for jstre = 1:nstre
                    for mstre = 1:nstre
                        for igaus = 1:ngaus
                            Smi = squeeze(tstres(mstre,igaus,istre,:));
                            Smj = squeeze(tstres(mstre,igaus,jstre,:));
                            factor = obj.computeVoigtFactor(mstre,nstre);
                            val = 1/V*(Smi.*Smj)'*dV(:,igaus)*factor;
                            Ch2(istre,jstre) =  Ch2(istre,jstre) + val;
                        end
                    end
                end
            end
            P = Shomog'*Ch2*Shomog;
        end
        
    end
end