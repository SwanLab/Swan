classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties (GetAccess = {?Physical_Problem, ?Element_Elastic, ?Element_Elastic_2D, ?Element_Elastic_3d, ?Element_Hyperelastic, ?Element_Elastic_Micro}, SetAccess = protected)
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
    end
    
    properties (GetAccess = {?Element_Elastic, ?Element_Elastic, ?Element_Thermal,?PhysicalVariables,?Element_Elastic_Micro}, SetAccess = {?Physical_Problem,?Element, ?Element_Elastic_Micro})
        B
    end
    
    methods (Access = ?Physical_Problem, Static)
        function element = create(ptype,pdim,dim,nelem,geometry,material,bc,dof)
            switch ptype
                case 'ELASTIC'
                    switch pdim
                        case '2D'
                            element = Element_Elastic_2D;
                            element.B = B2;
                        case '3D'
                            element = Element_Elastic_3D;
                            element.B = B3;
                    end
                case 'THERMAL'
                    element = Element_Thermal;
                    element.B = B_thermal;
                case 'HYPERELASTIC'
                    element = Element_Hyperelastic();
                otherwise
                    error('Invalid ptype.')
            end
            element.dim = dim;
            element.nunkn = dim.nunkn;
            element.nstre = dim.nstre;
            element.nelem = nelem;
            element.nnode = geometry.nnode;
            element.geometry = geometry;
            element.material = material;
            element.dof = dof;
            element.bc = bc;
            FextSupVol = element.computeExternalForces();
            element.assembleExternalForces(FextSupVol);
        end
    end
    
    
    methods (Access = {?Physical_Problem, ?Element})
        function FextSupVol = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext();
            FextVolumetric  = obj.computeVolumetricFext ();
            FextSupVol = FextSuperficial + FextVolumetric;
        end
        function obj = assembleExternalForces(obj,FextSupVol)
            % Assemble external forces
            obj.Fext = zeros(obj.dof.ndof,1);
            for i = 1:obj.nnode*obj.nunkn
                b = squeeze(FextSupVol(i,1,:));
                ind = obj.dof.idx(i,:);
                obj.Fext = obj.Fext + sparse(ind,1,b',obj.dof.ndof,1);
            end
            
            %Compute Global Puntual Forces (Not well-posed in FEM)
            if ~isempty(obj.bc.iN)
                FextPoint = zeros(obj.dof.ndof,1);
                FextPoint(obj.bc.iN) = obj.bc.neunodes(:,3);
                obj.Fext = obj.Fext + FextPoint;
            end
        end
        
        function b = AssembleVector(obj,b_elem)
            b = zeros(obj.dof.ndof,1);
            for i = 1:obj.nnode*obj.nunkn
                c = squeeze(b_elem(i,1,:));
                ind = obj.dof.idx(i,:);
                b = b + sparse(ind,1,c',obj.dof.ndof,1);
            end
        end
        
        
        
        function [A] = AssembleMatrix(obj,A_elem)
            % Assemble stiffness matrix
            A = sparse(obj.dof.ndof,obj.dof.ndof);
            for i = 1:obj.nnode*obj.nunkn
                for j = 1:obj.nnode*obj.nunkn
                    a = squeeze(A_elem(i,j,:));
                    A = A + sparse(obj.dof.idx(i,:),obj.dof.idx(j,:),a,obj.dof.ndof,obj.dof.ndof);
                end
            end
            A = 1/2 * (A + A');
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

