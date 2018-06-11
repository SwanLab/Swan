classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties (GetAccess = {?Physical_Problem, ?Element_Elastic, ?Element_Hyperelastic, ?Element_Elastic_2D, ?Element_Elastic_3d, ?Element_Hyperelastic, ?Element_Elastic_Micro}, SetAccess = {?Element_Hyperelastic, ?Element_Elastic, ?Element_Elastic_Micro, ?Element_Thermal, ?Physical_Problem})
        Fext
        nunkn
        nstre
        nelem
        nnode
        geometry
        material
        dof
        bc
        dim
        pdim
        coord
        fincr
        nincr
        cload
        uD
    end
    
    
    methods (Access = ?Physical_Problem, Static)
        function element = create(mesh,geometry,material,bc,dof,dim)
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
                            switch pdim
                                case {'2D','3D'}
%                                     element = Element_NeoHookean(geometry);
                                    element = Element_SaintVenant(geometry);
                            end
                        otherwise
                            error('Invalid ptype.')
                    end 
            end
            
            element.dim         = dim;
            element.nunkn       = dim.nunkn;
            element.nstre       = dim.nstre;
            element.nelem       = nelem;
            element.nnode       = geometry.nnode;
            element.geometry    = geometry;
            element.material    = material;
            element.dof         = dof;
            element.bc          = bc;
            element.coord       = mesh.coord;
            element.pdim        = pdim;
            element.assign_dirichlet_values()
            
            element.computeExternalForces();
            
            % Create force increment.
            element.Fext = 0.1*element.Fext;
%             element.cload = element.Fext;

            global test
            if exist('test','var') == 1
%                 fprintf('FEM-MAT-OO\n\nGeneral\n- force\t%f\n- incrm\t%d\n- model\t%s\n- elements\t%d\n\nMaterial\n- \x03bc\t%f\n- \x03bb\t%f\n\n',element.fincr(element.dof.neumann(1))*element.nincr,element.nincr,class(element),element.nelem,element.material.mu,element.material.lambda);
                fprintf('FEM-MAT-OO\n\nGeneral\n- force\t%f\n- incrm\t%d\n- model\t%s\n- elements\t%d\n\nMaterial\n- \x03bc\t%f\n- \x03bb\t%f\n\n',element.Fext(element.dof.neumann(1)),element.nincr,class(element),element.nelem,element.material.mu,element.material.lambda);
            end
        end
    end
    
    
    methods (Access = {?Physical_Problem, ?Element})
        function obj = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext();
            FextVolumetric  = obj.computeVolumetricFext ();
            FextSupVol = FextSuperficial + FextVolumetric;
            FextSupVol = obj.AssembleVector(FextSupVol);
            FextPoint = obj.computePunctualFext();
            obj.Fext = FextSupVol +  FextPoint;
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
        function b = AssembleVector(obj,b_elem)
            b = zeros(obj.dof.ndof,1);
            for i = 1:obj.nnode*obj.nunkn
                c = squeeze(b_elem(i,1,:));
                idof_elem = obj.dof.in_elem(i,:);
                b = b + sparse(idof_elem,1,c',obj.dof.ndof,1);
            end
        end
        
        
        % Matrix function
        function [A] = AssembleMatrix(obj,A_elem)
            A = sparse(obj.dof.ndof,obj.dof.ndof);
            for i = 1:obj.nnode*obj.nunkn
                for j = 1:obj.nnode*obj.nunkn
                    a = squeeze(A_elem(i,j,:));
                    A = A + sparse(obj.dof.in_elem(i,:),obj.dof.in_elem(j,:),a,obj.dof.ndof,obj.dof.ndof);
                end
            end
            A = 1/2 * (A + A');
        end
        
                
        function assign_dirichlet_values(obj)
            if ~isempty(obj.dof.dirichlet)
                obj.uD = obj.dof.dirichlet_values;
            else
                obj.uD = [];
            end
            
        end
        
        function R = compute_imposed_displacemet_force(obj,K)
            % Forces coming from imposed displacement
            if ~isempty(obj.dof.dirichlet)
                R = -K(:,obj.dof.dirichlet)*obj.uD;
            else
                R = zeros(obj.dof.ndof,1);
            end
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

