classdef Element_Elastic_Micro < handle
    
    properties (Access = protected)
        vstrain
    end
    
    methods (Access = public)
        
        function setVstrain(obj,s)
            obj.vstrain = s;
        end
        
    end
    
    methods (Access = protected)
        
        function variables = computeStressStrainAndCh(obj,variables)
            
            variables.stress_fluct = variables.stress;
            variables.strain_fluct = variables.strain;
            Cmat = obj.material.C;
            
            variables.stress = zeros(obj.quadrature.ngaus,obj.nstre,obj.nelem);
            variables.strain = zeros(obj.quadrature.ngaus,obj.nstre,obj.nelem);
            variables.stress_homog = zeros(obj.nstre,1);
            vol_dom = sum(sum(obj.geometry.dvolu));
            
            for igaus = 1:obj.quadrature.ngaus
                variables.strain(igaus,1:obj.nstre,:) = obj.vstrain.*ones(1,obj.nstre,obj.nelem) + variables.strain_fluct(igaus,1:obj.nstre,:);
                for istre = 1:obj.nstre
                    for jstre = 1:obj.nstre
                        variables.stress(igaus,istre,:) = squeeze(variables.stress(igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Cmat(istre,jstre,:))).* squeeze(variables.strain(igaus,jstre,:));
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
            eforce = zeros(obj.dof.nunkn*obj.nnode,1,obj.nelem);
            sigma = zeros(obj.nstre,1,obj.nelem);
            for igaus = 1:obj.quadrature.ngaus
                %                 Bmat = obj.computeB(obj.dof.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                Bmat = obj.computeB(igaus);
                for istre = 1:obj.nstre
                    for jstre = 1:obj.nstre
                        sigma(istre,:) = sigma(istre,:) + squeeze(Cmat(istre,jstre,:)*vstrain(jstre))';
                    end
                end
                for iv = 1:obj.nnode*obj.dof.nunkn
                    for istre = 1:obj.nstre
                        eforce(iv,:) = eforce(iv,:)+(squeeze(Bmat(istre,iv,:)).*sigma(istre,:)'.*obj.geometry.dvolu(:,igaus))';
                    end
                end
            end
            F = -eforce;
        end
        

 
    end
end