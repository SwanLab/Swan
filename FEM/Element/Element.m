classdef Element < handle
    %Element Summary of this class goes here
    
    properties %(GetAccess = {?Physical_Problem, ?Element_Elastic, ?Element_Thermal, ?Element_Hyperelastic, ?Element_Elastic_2D, ?Element_Elastic_3d, ?Element_Hyperelastic, ?Element_Elastic_Micro}, SetAccess = protected)
        nstre
        nelem
        nnode
        geometry
        quadrature
        material
        dof
        uD
        nfields
    end
    
%     properties (Abstract, Access = protected)
%         scale
%     end
    
    properties (Access = protected)
        bcApplier        
    end

    methods (Static)        
        function obj = Element(geometry,material,dof,scale)            
            obj.nelem = geometry(1).interpolation.nelem;
            obj.nfields = geometry.nfields;
            for ifield=1:obj.nfields
                obj.nnode(ifield) = geometry(ifield).interpolation.nnode;
            end
            obj.geometry = geometry;
            obj.quadrature = Quadrature.set(geometry(1).type);
            obj.material = material;
            obj.dof = dof;
            obj.bcApplier = BoundaryConditionsApplier.create(obj.nfields,obj.dof,scale);
            obj.assign_dirichlet_values();
        end
    end
    
    methods
        
        function bc = getBcApplier(obj)
            bc = obj.bcApplier;            
        end
        
        function [r,dr] = computeResidual(obj,x)
            % !! Currently unused !!
            % *************************************************************
            % Compute
            % - residual: r = LHS*x - RHS
            % - residual derivative: dr = LHS
            % *************************************************************
            
            RHS = obj.computeRHS;
            LHS = obj.computeLHS*x;
            
            r = LHS*x - RHS;
            dr = LHS;
        end
        
        function Fext = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext;
            FextVolumetric  = obj.computeVolumetricFext;
            FextSupVol = {FextSuperficial + FextVolumetric};
            FextSupVol = obj.AssembleVector(FextSupVol);
            FextPoint = obj.computePunctualFext;
            Fext = FextSupVol +  FextPoint;
        end
        
        % *****************************************************************
        % Assembling Functions
        %******************************************************************
        function FextPoint = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            FextPoint = zeros(obj.dof.ndof,1);
            if ~isempty(obj.dof.neumann)
                FextPoint(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end
        
        % Vector function
        function b = AssembleVector(obj,b_elem_cell)
            for ifield = 1:obj.nfields
                b_elem = b_elem_cell{ifield,1};
                b = zeros(obj.dof.ndof(ifield),1);
                for i = 1:obj.geometry(ifield).interpolation.nnode*obj.dof.nunkn(ifield)
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
            nunkn1 = obj.dof.nunkn(ifield);
            nnode1 = obj.geometry(ifield).interpolation.nnode;
            nunkn2 = obj.dof.nunkn(jfield);
            nnode2 = obj.geometry(jfield).interpolation.nnode;
            col = obj.dof.ndof(jfield);
            row = obj.dof.ndof(ifield);
        end
        
        function assign_dirichlet_values(obj)
            for ifield = 1:obj.nfields
                if ~isempty(obj.dof.dirichlet{ifield})
                    obj.uD{ifield} = obj.dof.dirichlet_values{ifield};
                else
                    obj.uD = {[]};
                end
            end
        end
        
        function R = compute_imposed_displacement_force(obj,K)
            % Forces coming from imposed displacement
            dirApplier = DirichletConditionsApplier(obj.nfields,obj.dof);
            [dirichlet,uD,~] = dirApplier.compute_global_dirichlet_free_uD();
            if ~isempty(dirichlet)
                R = -K(:,dirichlet)*uD;
            else
                R = zeros(sum(obj.dof.ndof(:)),1);
            end
        end
        
 
    end
    
    methods (Access = protected)
       
        function init(obj,geo)
            
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

