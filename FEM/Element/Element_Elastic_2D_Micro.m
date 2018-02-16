classdef Element_Elastic_2D_Micro < Element_Elastic_2D
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
        vstrain
    end
    
    methods 
        
        
        function variables = computeVars(obj,uL)
            variables = computeVars@Element_Elastic_2D(obj,uL);
            
            variables.stress_fluct = variables.stress;
            variables.strain_fluct = variables.strain;
            Cmat = obj.material.C;
            
            variables.stress = zeros(obj.geometry.ngaus,obj.nstre,obj.nelem);
            variables.strain = zeros(obj.geometry.ngaus,obj.nstre,obj.nelem);
            variables.stress_homog = zeros(obj.dim.nstre,1);
            vol_dom = sum(sum(obj.geometry.dvolu));
            
            for igaus=1:obj.geometry.ngaus
                variables.strain(igaus,1:obj.nstre,:) = obj.vstrain.*ones(1,obj.nstre,obj.nelem) + variables.strain_fluct(igaus,1:obj.nstre,:);
                for istre=1:obj.dim.nstre
                    for jstre=1:obj.dim.nstre
                        variables.stress(igaus,istre,:) = squeeze(variables.stress(igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Cmat(istre,jstre,:))).* squeeze(variables.strain(igaus,jstre,:));
                    end
                end
                % contribucion a la C homogeneizada
                for istre=1:obj.nstre
                    variables.stress_homog(istre) = variables.stress_homog(istre) +  1/vol_dom *(squeeze(variables.stress(igaus,istre,:)))'*obj.geometry.dvolu(:,igaus);
                end
            end
        end
        
        end
    
   
    
    methods (Access = {?Physical_Problem, ?Element})
%         function obj = Element_Elastic_Micro
%             obj.B = B2;
%         end
%         
%         function obj = computeRHS(obj,vstrain)
%             computeRHS@Element(obj);
%             RHSStrain = obj.computeStrainRHS(vstrain);
%             obj.RHS = obj.RHS + RHSStrain;             
%         end
    end
    
    methods   %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?Element})
        function [r,dr] = computeResidual(obj,x)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************

            
            % Compute stiffness matrix
            [K] = obj.computeStiffnessMatrix();
            
            % Assemble
            [K] = obj.AssembleMatrix(K);
           
            %Set fext
            Fext = obj.computeExternalForces();
            R = obj.compute_imposed_displacemet_force(K); 
            obj.fext = Fext+R;
            
            
            vF = obj.dof.vF;
            vP = obj.dof.vP;
            vQ = obj.dof.vQ;
            vI = setdiff(vF,vP);
                  
            K_II = K(vI,vI);
            K_IP = K(vI,vP) + K(vI,vQ); %Grouping P and Q nodal values
            K_PI = K(vP,vI) + K(vQ,vI); % Adding P  and Q equation 
            K_PP = K(vP,vP) + K(vP,vQ) + K(vQ,vP) + K(vQ,vQ); % Adding and grouping
            
            Kred = [K_II K_IP; K_PI K_PP];
            
            fext_I = obj.fext(vI);
            fext_P = obj.fext(vP)+obj.fext(vQ);
            
            fext_red = [fext_I; fext_P];
            
            fint_red = Kred*x;
            
            r = fint_red - fext_red;
            dr = Kred;

            
 
        end
           
        
    end
    
    methods (Access = protected)
        function u = compute_displacements(obj,usol)
            u = zeros(obj.dof.ndof,1);
            u(obj.dof.vF) = usol;
            u(obj.dof.vQ) = u(obj.dof.vP);
            u(obj.dof.vD) = obj.uD;
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = computeVolumetricFext@Element_Elastic(obj,bc);
            F_def = obj.computeStrainRHS(obj.vstrain);
            FextVolumetric = FextVolumetric + F_def;
        end

    end
    
    
    methods (Access = private)
        function F = computeStrainRHS(obj,vstrain)
            Cmat = obj.material.C;            
            eforce = zeros(obj.nunkn*obj.nnode,1,obj.nelem);
            sigma=zeros(obj.nstre,1,obj.nelem);
            for igaus=1:obj.geometry.ngaus
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                for istre=1:obj.nstre
                    for jstre=1:obj.nstre
                        sigma(istre,:) = sigma(istre,:) + squeeze(Cmat(istre,jstre,:)*vstrain(jstre))';
                    end
                end                
                for iv=1:obj.nnode*obj.nunkn
                    for istre=1:obj.nstre
                        eforce(iv,:)=eforce(iv,:)+(squeeze(Bmat(istre,iv,:)).*sigma(istre,:)'.*obj.geometry.dvolu(:,igaus))';
                    end
                end
            end
            F = -eforce;
        end
    end
end
