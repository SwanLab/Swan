classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties %(GetAccess = {?Physical_Problem, ?Element_Elastic, ?Element_Thermal, ?Element_Hyperelastic, ?Element_Elastic_2D, ?Element_Elastic_3d, ?Element_Hyperelastic, ?Element_Elastic_Micro}, SetAccess = protected)
        nunkn
        nstre
        nelem
        nnode
        geometry
        material
        dof
        bc
        dim
        uD
        nfields
    end
    
    
    methods (Access = ?Physical_Problem, Static)
        function element = create(mesh,geometry,material,bc,dof,dim,nfields)
            
            nelem = mesh.nelem;
            ptype = mesh.ptype;
            pdim = mesh.pdim;
            
            switch mesh.scale
                
                case 'MICRO'
                    element = Element_Elastic_2D_Micro;
                case 'MACRO'
                    switch ptype
                        case 'ELASTIC'
                            switch pdim
                                case '2D'
                                    element = Element_Elastic_2D;
                                case '3D'
                                    element = Element_Elastic_3D;
                            end
                        case 'THERMAL'
                            element = Element_Thermal;
                        case 'HYPERELASTIC'
                            element = Element_Hyperelastic();
                        case 'Stokes'
                            element = Element_Stokes;
                        otherwise
                            error('Invalid ptype.')
                    end 
            end
            
            element.nfields = nfields;
            for ifield=1:nfields
                element.nunkn(ifield) = dim.nunkn(ifield);
                element.nnode(ifield) = geometry(ifield).nnode;
            end
            element.dim = dim;
            element.nstre = dim.nstre;
            element.nelem = nelem;
            element.geometry = geometry;
            element.material = material;
            element.dof = dof;
            element.bc = bc;
            element.assign_dirichlet_values()
        end
    end
    
    
    methods (Access = {?Physical_Problem, ?Element})
        function Fext = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext();
            FextVolumetric  = obj.computeVolumetricFext ();
            FextSupVol = {FextSuperficial + FextVolumetric};
            FextSupVol = obj.AssembleVector(FextSupVol);
            FextPoint = obj.computePunctualFext();
            Fext = FextSupVol +  FextPoint;
        end
        
        % *****************************************************************
        % Assembling Functions
        %******************************************************************
         function FextPoint = computePunctualFext(obj)    
            %Compute Global Puntual Forces (Not well-posed in FEM)
            if ~isempty(obj.dof.neumann)
                FextPoint = zeros(obj.dof.ndof,1);
                FextPoint(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end
        
        % Vector function
        function b = AssembleVector(obj,b_elem_cell)
           for ifield = 1:obj.nfields
                b_elem = b_elem_cell{ifield,1};
                b = zeros(obj.dof.ndof(ifield),1);
                for i = 1:obj.nnode(ifield)*obj.nunkn(ifield)
                    c = squeeze(b_elem(i,1,:));
                    idof_elem = obj.dof.in_elem{ifield}(i,:);
                    b = b + sparse(idof_elem,1,c',obj.dof.ndof(ifield),1);
                end
                b_global{ifield,1} = b; 
           end
           b=cell2mat(b_global);
        end
        
        
        % Matrix function
        function [A] = AssembleMatrix(obj,A_elem,ifield,jfield)
            
              [idx1,idx2,nunkn1,nunkn2,nnode1,nnode2,col,row] = obj.get_assemble_parameters(ifield,jfield);
                    
                    A = sparse(row,col);
                    for i = 1:nnode1*nunkn1
                        for j = 1:nnode2*nunkn2
                            a = squeeze(A_elem(i,j,:));
                            A = A + sparse(idx1(i,:),idx2(j,:),a,row,col);
                        end
                    end
                   
                    if ifield == 1 && jfield == 1
                         A = 1/2 * (A + A');
                    end
        end
        
        function [idx1,idx2,nunkn1,nunkn2,nnode1,nnode2,col,row] = get_assemble_parameters(obj,ifield,jfield)
                    idx1 = obj.dof.in_elem{ifield};
                    idx2 = obj.dof.in_elem{jfield};
                    nunkn1 = obj.dim.nunkn(ifield);
                    nnode1 = obj.geometry(ifield).nnode;
                    nunkn2 = obj.dim.nunkn(jfield);
                    nnode2 = obj.geometry(jfield).nnode;
                    col = obj.dof.ndof(jfield);
                    row = obj.dof.ndof(ifield);
        end
                
        function assign_dirichlet_values(obj)
            for ifield = 1:obj.nfields
                if ~isempty(obj.dof.dirichlet{ifield})
                    obj.uD{ifield} = obj.dof.dirichlet_values{ifield};
                else
                    obj.uD = [];
                end
            end
  
        end
        
        function R = compute_imposed_displacemet_force(obj,K)
            % Forces coming from imposed displacement
             [dirichlet,uD,~] = obj.compute_global_dirichlet_free_uD;
                if ~isempty(dirichlet)
                    R = -K(:,dirichlet)*uD;
                else
                    R = zeros(sum(obj.dof.ndof(:)),1);
                end
            

        end
        
        function [dirichlet,uD,free] = compute_global_dirichlet_free_uD(obj)
            global_ndof=0;
            for ifield=1:obj.nfields
                dirichlet{ifield,1} = obj.dof.dirichlet{ifield}+global_ndof;
                uD{ifield,1} = obj.uD{ifield};
                free{ifield,1} = obj.dof.free{ifield}' + global_ndof;
                global_ndof=global_ndof+obj.dof.ndof(ifield);
            end
            uD = cell2mat(uD);
            dirichlet = cell2mat(dirichlet);
            free = cell2mat(free);
        end

         function Ared = full_matrix_2_reduced_matrix(obj,A)
             [~,~,free] = obj.compute_global_dirichlet_free_uD;
 
                Ared = A(free,free);

        end
        
        function b_red = full_vector_2_reduced_vector(obj,b)
            [~,~,free] = obj.compute_global_dirichlet_free_uD;
           
                b_red = b(free);
        
        end
        
        function b = reduced_vector_2_full_vector(obj,bfree)
           [dirichlet,uD,free] = obj.compute_global_dirichlet_free_uD; 
           ndof = sum(obj.dof.ndof);
            
                b = zeros(ndof,1);
                b(free) = bfree;
                b(dirichlet) = uD;
       
        end
        
    end
    
    methods (Abstract, Access = protected)
        r = computeSuperficialFext(obj)
        r = computeVolumetricFext(obj)
    end
    
%     methods (Abstract, Access = {protected, ?Physical_Problem})
%             [r,dr] = computeResidual(x)
%     end
end

