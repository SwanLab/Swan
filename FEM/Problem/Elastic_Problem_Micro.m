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
            for istre=1:nstre
                obj.element.setVstrain(vstrain(istre,:));
                obj.computeVariables;
                Chomog(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstress(istre,:,:,:) = obj.variables.stress;
            end
            obj.variables.Chomog = Chomog;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
        end
        
        function P2 = computeAmplificator(obj)
            obj.variables.Ptensor = rand(3,3);
            [Chomog,tstrain,~] = obj.computeChomog();
             Ct = SymmetricFourthOrderPlaneStressVoigtTensor();
             Ct.setValue(Chomog);
             St = Inverter.invert(Ct);
             Shomog = St.getValue();

            identStrain = ones(size(tstrain));
            corrTens = identStrain + tstrain;
            W = corrTens;
            
            V = sum(sum(obj.geometry.dvolu));
            ngaus = obj.element.quadrature.ngaus;
            dV    = obj.element.geometry.dvolu;
            C     = obj.element.material.C;
            nstre = obj.element.getNstre();
            
            P2 = zeros(size(Chomog));

            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        for mstre = 1:nstre
                            for kstre = 1:nstre
                                for lstre = 1:nstre
                                    pim = C(istre,kstre,igaus,:).*W(kstre,igaus,lstre,:)*Shomog(lstre,mstre);
                                    Pim = 1/V*squeeze(pim)'*dV(:,igaus);
                                    pmj = C(mstre,kstre,igaus,:).*W(kstre,igaus,lstre,:)*Shomog(lstre,jstre);
                                    Pmj = 1/V*squeeze(pmj)'*dV(:,igaus);
                                    P2(istre,jstre) = P2(istre,jstre) + Pim*Pmj;
                                end
                            end
                        end
                    end
                end
            end
            
            
        end
        
    end
end