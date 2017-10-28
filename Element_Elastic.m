classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = Element_Elastic(ndim)
            switch ndim
                case 2
                    % !! B should be property of Element_Elastic or Element !!
                    obj.B = B2;
                case 3
                    obj.B = B3;
            end
        end
        
        function obj = computeLHS(obj,nunkn,nelem,geometry,material)
            
            Ke = zeros(nunkn*geometry.nnode,nunkn*geometry.nnode,nelem);
            % Elastic matrix
            Cmat = material.C;
            
            for igauss=1:geometry.ngaus
                % Strain-displacement matrix
                [obj.B, Bmat] = obj.B.computeB(nunkn,nelem,geometry.nnode,geometry.cartDeriv(:,:,:,igauss));
                
                % Compute Ke
                for i = 1:nelem
                    Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*Cmat(:,:,i)*...
                        Bmat(:,:,i)*geometry.area(i,igauss);
                end
            end
            obj.LHS = Ke;
        end
        
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)
            Fext = zeros(nnode*nunkn,1,nelem);
            for i = 1:length(bc.iN)
                for j = 1:nelem
                    ind = find(idx(:,j) == bc.iN(i));
                    if ~isempty(ind)
                        Fext(ind,:,j) = bc.neunodes(i,3);
                    end
                    % clear ind
                    ind = [];
                end
            end
            obj.RHS = Fext;
        end
        
    end
    
end
