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
        elem_obj
        row_index
        column_index 
        Cmat
        Bmat
        CtimesB
        inodes
        icomps
        InitialNonAssembledIndex
        FinalNonAssembledIndex
        Ke
        connectivities
        Subintegrated_StifMat
        StifMat
        dvolum
        ngaus
        numberOfGlobalDofs
        rowGlobalDofs
        columnGlobalDofs
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
            obj.Cmat = elem_obj.material.C;
            obj.elem_obj = elem_obj;
            obj.inodes=reshape(repmat(1:obj.nnode,obj.nunkn,1),1,[]); 
            obj.icomps=repmat(1:obj.nunkn,1,obj.nnode);
        end
        
        function generate(obj)
            
            obj.initilize_StifMat()
            obj.compute_dvolum() 
                           
            for igaus=1:obj.ngaus
                obj.initialize_row_and_column_index();
                obj.computeBmatrix(igaus);
                obj.compute_CtimesBmatrix();
                obj.computeNonDiagonalEntries(igaus);            
                obj.computeDiagonalEntries(igaus);
                obj.assemble_matrix()
                obj.add_matrix()
            end
            obj.symmetrizeStiffMat()
        end
        
        function K = getStiffMatrix(obj)
            K = obj.StifMat;
        end
        
    end
    
    methods (Access = private)
        
        function initilize_StifMat(obj)
            obj.StifMat = sparse(obj.number_of_global_dofs,obj.number_of_global_dofs);
        end
        
        
        function compute_dvolum(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            obj.ngaus = obj.quadrature.ngaus;
            obj.dvolum = obj.geometry.dvolu;
        end
        
        function computeBmatrix(obj,igaus)
            obj.Bmat = obj.elem_obj.computeB(igaus);
        end
        
        function compute_CtimesBmatrix(obj)
            CB = zeros(obj.nstre,obj.number_of_elemental_dofs,obj.nelem);
            for i=1:obj.nstre
                CB(i,:,:) = sum(repmat(permute(obj.Cmat(i,:,:),[2,1,3]),1,obj.number_of_elemental_dofs,1) .* obj.Bmat,1);
            end
            obj.CtimesB = CB;
        end
        
        function initialize_row_and_column_index(obj)
            obj.row_index = obj.initialize_index();
            obj.column_index = obj.initialize_index();
        end
        
        
        function index = initialize_index(obj)
            index = zeros(obj.number_of_elemental_dofs*obj.number_of_elemental_dofs*obj.nelem,1);
        end
        
        function  computeDiagonalEntries(obj,igaus)
                for idof=1:obj.number_of_elemental_dofs
                    obj.obtainGlobalDofsForSpecificRowEntry(idof);
                    obj.obtainNumberOfGlobalDofs();
                    obj.FinalNonAssembledIndex = obj.InitialNonAssembledIndex + obj.numberOfGlobalDofs -1;
                        
                    k_ij = obj.computeStiffEntry(idof,idof,igaus);
                    
                    index3 = obj.compute_index3();
                    
                    obj.row_index(index3,1)    =  obj.rowGlobalDofs ;
                    obj.column_index(index3,1) =  obj.rowGlobalDofs ;
                    obj.Ke(index3,1) =  k_ij;
                    obj.InitialNonAssembledIndex = obj.FinalNonAssembledIndex + 1;
                end
        end
        
        
        function computeNonDiagonalEntries(obj,igaus)
            obj.InitialNonAssembledIndex=1;   
            obj.initialize_Ke()
            for idof=1:obj.number_of_elemental_dofs
                    obj.obtainGlobalDofsForSpecificRowEntry(idof);
                    obj.obtainNumberOfGlobalDofs();
                    
                    for jdof=1:idof-1
                        obj.FinalNonAssembledIndex = obj.InitialNonAssembledIndex + 2*obj.numberOfGlobalDofs -1;
                        obj.obtainGlobalDofsForSpecificColumnEntry(jdof)                        
                        k_ij = obj.computeStiffEntry(idof,jdof,igaus);
                        
                        index2 = obj.compute_index2();
                        
                        obj.row_index(index2,1) =    [ obj.rowGlobalDofs    ; obj.columnGlobalDofs];
                        obj.column_index(index2,1) = [ obj.columnGlobalDofs ; obj.rowGlobalDofs];
                       
                        obj.Ke(index2,1) = [k_ij ; k_ij ];
                        
                        obj.InitialNonAssembledIndex = obj.FinalNonAssembledIndex + 1;
                    end
            end    
        end
        
        function initialize_Ke(obj)
            obj.Ke = zeros(obj.number_of_elemental_dofs*obj.number_of_elemental_dofs*obj.nelem,1);
        end
        
        
        function K = computeStiffEntry(obj,idof,jdof,igaus)
            BCB = squeeze(sum(obj.Bmat(:,idof,:) .* obj.CtimesB(:,jdof,:),1));
            K = obj.dvolum(:,igaus).*BCB;
        end
        
        function obtainGlobalDofsForSpecificRowEntry(obj,idof)
            obj.rowGlobalDofs = obj.compute_global_dofs(idof);
        end
        
        function obtainGlobalDofsForSpecificColumnEntry(obj,idof)
            obj.columnGlobalDofs = obj.compute_global_dofs(idof);
        end
        
        function obtainNumberOfGlobalDofs(obj)
            obj.numberOfGlobalDofs = length(obj.rowGlobalDofs);

        end
        
        function  global_dofs = compute_global_dofs(obj,idof)
            global_dofs = obj.nunkn*(obj.connectivities(:,obj.inodes(idof))-1)+obj.icomps(idof);
        end
        
        function index2 = compute_index2(obj)
            index2 = obj.InitialNonAssembledIndex:obj.FinalNonAssembledIndex;
        end

        function index3 = compute_index3(obj)
             index3 = obj.InitialNonAssembledIndex:obj.FinalNonAssembledIndex;
        end
        
        function assemble_matrix(obj)
             obj.Subintegrated_StifMat = sparse(obj.row_index,obj.column_index,obj.Ke,obj.number_of_global_dofs,obj.number_of_global_dofs);
        end
        
        function add_matrix(obj)
            obj.StifMat = obj.StifMat + obj.Subintegrated_StifMat;
        end
        
        function symmetrizeStiffMat(obj)
            obj.StifMat = 1/2 * (obj.StifMat + obj.StifMat');
        end
        

    end
    
end

