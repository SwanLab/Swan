classdef StifnessMatrxiGenerator < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    properties (Access = private)
        nnode
        nunkn
        nstre
        quadrature
        interpolation_u
        geometry
        number_of_global_dofs
        number_of_elemental_dofs
        nelem
        material
        elem_obj
        row_index
        column_index 
        Cmat
        Bmat
        CtimesB
        inodes
        icomps
        someindex
        Ke
        connectivities
    end
    
    methods
        
        function obj= StifnessMatrxiGenerator(elem_obj)
            obj.nnode = elem_obj.nnode;
            obj.nunkn = elem_obj.dof.nunkn;
            obj.nstre = elem_obj.nstre;
            obj.number_of_global_dofs  = elem_obj.dof.ndof;
            obj.number_of_elemental_dofs = obj.nnode*obj.nunkn;
            obj.nnode = elem_obj.nnode;
            obj.nelem = elem_obj.nelem;
            
            obj.quadrature = elem_obj.quadrature;
            obj.interpolation_u = elem_obj.interpolation_u;
            obj.geometry = elem_obj.geometry;
            obj.connectivities = elem_obj.geometry.interpolation.T;
            obj.material = elem_obj.material;
            obj.elem_obj = elem_obj;
        end
        
        function K = generate(obj)
            
            obj.quadrature.computeQuadrature('LINEAR');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            % Stiffness matrix
            StifMat = sparse(obj.number_of_global_dofs,obj.number_of_global_dofs);
            obj.inodes=reshape(repmat(1:obj.nnode,obj.nunkn,1),1,[]); 
            obj.icomps=repmat(1:obj.nunkn,1,obj.nnode);
                
                
            % Elastic matrix
            obj.Cmat = obj.material.C;
            for igaus=1:obj.quadrature.ngaus
                
                obj.Bmat = obj.elem_obj.computeB(igaus);
                obj.compute_CtimesBmatrix();
                obj.initialize_row_and_column_index_to_zero();
                
                obj.computeNonDiagonalEntries(igaus);            
                obj.computeDiagonalEntries(igaus);
                
                StifMat_of_iguass = compute_assembled_matrix();
                StifMat = StifMat + StifMat_of_iguass;
            end
            K = 1/2 * (StifMat + StifMat');           
            
        end
    end
    
    methods (Access = private)
        
            function compute_CtimesBmatrix(obj)
            CB =zeros(obj.nstre,obj.number_of_elemental_dofs,obj.nelem);
            for i=1:obj.nstre
                CB(i,:,:) = sum(repmat(permute(obj.Cmat(i,:,:),[2,1,3]),1,obj.number_of_elemental_dofs,1) .* obj.Bmat,1);
            end
            obj.CtimesB = CB;
        end
        
        function initialize_row_and_column_index_to_zero(obj)
            obj.row_index = obj.initialize_index();
            obj.column_index = obj.initialize_index();
        end
        
        function index = initialize_index(obj)
            index = zeros(obj.number_of_elemental_dofs*obj.number_of_elemental_dofs*obj.nelem,1);
        end
        
        function  computeDiagonalEntries(obj,igaus)
                for idof=1:obj.number_of_elemental_dofs
                    it = obj.compute_global_dofs(idof);
                    
                    k_ij=squeeze(sum(obj.Bmat(:,idof,:) .* obj.CtimesB(:,idof,:),1));
                    
                    index3 = obj.compute_index3(it);
                    
                    obj.row_index(index3,1) =  it ;
                    obj.column_index(index3,1) =  it ;
                    obj.Ke(index3,1) =  obj.geometry.dvolu(:,igaus).*k_ij ;
                    obj.someindex = obj.someindex + length(it) ;
                end
        end
        
        
        function computeNonDiagonalEntries(obj,igaus)
            obj.someindex=1;   
            obj.Ke = zeros(obj.number_of_elemental_dofs*obj.number_of_elemental_dofs*obj.nelem,1);
            for idof=1:obj.number_of_elemental_dofs
                    it= obj.compute_global_dofs(idof);
                    for jdof=1:idof-1
                        jt = obj.compute_global_dofs(jdof);
                        
                        k_ij=squeeze(sum(obj.Bmat(:,idof,:) .* obj.CtimesB(:,jdof,:),1));
                        
                        index2 = obj.compute_index2(it);
                        
                        obj.row_index(index2,1) = [ it ; jt];
                        obj.column_index(index2,1) = [ jt ; it];
                        obj.Ke(index2,1) = [obj.geometry.dvolu(:,igaus).*k_ij ; obj.geometry.dvolu(:,igaus).*k_ij ];
                        obj.someindex = obj.someindex + 2*length(it);
                    end
            end    
        end
        
        
        function  global_dofs = compute_global_dofs(obj,idof)
            global_dofs = obj.nunkn*(obj.connectivities(:,obj.inodes(idof))-1)+obj.icomps(idof);
        end
        
        function index2 = compute_index2(obj,it)
            index2 = obj.someindex:obj.someindex+2*length(it)-1;
        end

        function index3 = compute_index3(obj,it)
             index3 = obj.someindex:obj.someindex+length(it)-1;
        end
        
        function Assembled_matrix = compute_assembled_matrix(obj)
            Assembled_matrix = sparse(obj.row_index,obj.column_index,obj.Ke,obj.number_of_global_dofs,obj.number_of_global_dofs);
        end
    end
    
end

