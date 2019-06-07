classdef Element_DiffReact_Micro < Element_DiffReact
    
    methods
        
        function FextPoint = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-fashioned in FEM)
            if ~isempty(obj.dof.neumann)
                FextPoint = obj.reducedToFullVector(obj.dof.neumann_values);
            end
        end
        
        function Ared = fullToReducedMatrix(obj,A)                
            vF = obj.dof.free;
            vP = obj.dof.periodic_free;
            vQ = obj.dof.periodic_constrained;
            vI = setdiff(vF{1},vP);
            
            A_II = A(vI,vI);
            A_IP = A(vI,vP) + A(vI,vQ); %Grouping P and Q nodal values
            A_PI = A(vP,vI) + A(vQ,vI); % Adding P  and Q equation
            A_PP = A(vP,vP) + A(vP,vQ) + A(vQ,vP) + A(vQ,vQ); % Adding and grouping
            
            Ared = [A_II, A_IP; A_PI, A_PP];
        end
        
        function b_red = fullToReducedVector(obj,b)
            vF = obj.dof.free{1};
            vP = obj.dof.periodic_free;
            vQ = obj.dof.periodic_constrained;
            vI = setdiff(vF,vP);
            
            b_I = b(vI);
            b_P = b(vP) + b(vQ);
            b_red = [b_I; b_P];
        end
        
        function b = reducedToFullVector(obj,bfree)
            % HEAD
            % b = zeros(obj.dof.ndof,1);
            % b(obj.dof.free{1}) = bfree;
            % b(obj.dof.dirichlet{1}) = obj.dof.dirichlet_values{1};
            
            % MASTER
            
            vF = obj.dof.free;
            vP = obj.dof.periodic_free;
            vI = setdiff(vF{1},vP);
            
            b = zeros(obj.dof.ndof,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(obj.dof.periodic_free) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(obj.dof.periodic_constrained) = b(obj.dof.periodic_free);
            
            %             b(obj.dof.free) = bfree;
            %             b(obj.dof.dirichlet) = obj.dof.dirichlet_values;
            %             b(obj.dof.periodic_constrained) = b(obj.dof.periodic_free);
        end
    end
    
    methods (Access = private)
        function F = computeStrainRHS(obj,vstrain)            
            Cmat = obj.material.C;
            eforce = zeros(obj.dof.nunkn*obj.nnode,1,obj.nelem);
            sigma = zeros(obj.nstre,1,obj.nelem);
            for igaus = 1:obj.quadrature.ngaus
                Bmat = obj.computeB(obj.dof.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                for istre = 1:obj.nstre
                    for jstre = 1:obj.nstre
                        sigma(istre,:) = sigma(istre,:) + squeeze(Cmat(istre,jstre)*vstrain(jstre))';
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
