classdef Element_Elastic_Micro < Element_Elastic
    
    properties (Access = protected)
        vstrain
    end
    
    methods (Access = public)
        
        function setVstrain(obj,s)
            obj.vstrain = s;
        end
        
        function variables = computeVars(obj,uL)
            variables = obj.computeVars@Element_Elastic(uL);
            variables = obj.computeStressStrainAndCh(variables);
        end
        
    end
    
    methods (Access = protected)
            
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = obj.computeVolumetricFext@Element_Elastic();
            F_def = obj.computeStrainRHS(obj.vstrain);
            FextVolumetric = FextVolumetric + F_def;
        end
        
    end
    
    methods (Access = private)
        
        function variables = computeStressStrainAndCh(obj,variables)
            
            variables.stress_fluct = variables.stress;
            variables.strain_fluct = variables.strain;
            Cmat = obj.material.C;
            
            variables.stress = zeros(obj.quadrature.ngaus,obj.nstre,obj.nelem);
            variables.strain = zeros(obj.quadrature.ngaus,obj.nstre,obj.nelem);
            variables.stress_homog = zeros(obj.nstre,1);
            vol_dom = 1;%sum(sum(obj.geometry.dvolu));
            
            for igaus = 1:obj.quadrature.ngaus
                variables.strain(igaus,1:obj.nstre,:) = obj.vstrain.*ones(1,obj.nstre,obj.nelem) + variables.strain_fluct(igaus,1:obj.nstre,:);
                for istre = 1:obj.nstre
                    for jstre = 1:obj.nstre
                        C  = squeeze(squeeze(Cmat(istre,jstre,:,igaus)));
                        variables.stress(igaus,istre,:) = squeeze(variables.stress(igaus,istre,:)) + 1/vol_dom*C.* squeeze(variables.strain(igaus,jstre,:));
                    end
                end
                % contribucion a la C homogeneizada
                for istre = 1:obj.nstre
                    variables.stress_homog(istre) = variables.stress_homog(istre) +  1/vol_dom *(squeeze(variables.stress(igaus,istre,:)))'*obj.geometry.dvolu(:,igaus);
                end
                
            end
            
        end
        
        function F = computeStrainRHS(obj,vstrain)
            Cmat = obj.material.C;
            ngaus = obj.quadrature.ngaus;
            eforce = zeros(obj.dof.nunkn*obj.nnode,ngaus,obj.nelem);
            sigma = zeros(obj.nstre,ngaus,obj.nelem);
            for igaus = 1:ngaus
                Bmat    = obj.computeB(igaus);
                dV(:,1) = obj.geometry.dvolu(:,igaus);                
                for istre = 1:obj.nstre
                    for jstre = 1:obj.nstre
                        Cij = squeeze(Cmat(istre,jstre,:,igaus));
                        vj  = vstrain(jstre);
                        si  = squeeze(sigma(istre,igaus,:));
                        sigma(istre,igaus,:) = si + Cij*vj;
                    end
                end
                for iv = 1:obj.nnode*obj.dof.nunkn
                    for istre = 1:obj.nstre
                        Biv_i = squeeze(Bmat(istre,iv,:));
                        si    = squeeze(sigma(istre,igaus,:));
                        Fiv   = squeeze(eforce(iv,igaus,:));
                        eforce(iv,igaus,:) = Fiv + Biv_i.*si.*dV;
                    end
                end
            end
            F = -eforce;
        end
        

 
    end
end