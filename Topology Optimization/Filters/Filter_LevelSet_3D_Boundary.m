classdef Filter_LevelSet_3D_Boundary < Filter_LevelSet_3D & Filter_LevelSet_Boundary
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function setInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Triangle_Linear(obj.unfitted_mesh);
        end
        
        function S = computeSurface(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            M2 = obj.computeRHS(x,ones(size(x)));
            S = sum(M2);
            
            filter2D = Filter_P1_LevelSet_2D_Interior;
            filter2D.loadProblem(obj.diffReacProb.problemID,'MACRO');
            
            for idime = 1:obj.mesh.ndim
                for iside = 1:2
                    face_mesh = obj.createFaceMesh(idime,iside);
                    obj.unfitted_mesh = Mesh_Unfitted_2D_Interior(face_mesh,obj.interpolation_unfitted);
                    obj.unfitted_mesh.computeMesh(x);
                    if ~isempty(obj.unfitted_mesh.fitted_cut_cells)
                        M2 = obj.computeRHS(x,ones(size(x)));
                        S = S + M2;
                    end
                end
            end
        end
    end
    
    methods ( Access = private)
        function face_mesh = createFaceMesh(obj,idime,iside)
            dimen = [1 2 3];
            face_connec = obj.getFaceConnectivities(idime,iside);
            face_mesh = Mesh;
            
            face_mesh = face_mesh.create(obj.mesh.coord(:,dimen(dimen ~= idime)),face_connec);
        end
        
        function face_connec = getFaceConnectivities(obj,idime,iside)
            if iside == 1
                L = min(obj.mesh.coord(:,idime));
            else
                L = max(obj.mesh.coord(:,idime));
            end
            
            valid_coord = obj.mesh.coord(:,idime) == L;
            coord_indexes = find(valid_coord);
            switch obj.nnode
                case 4
                    new_nnode = 3;
                case 8
                    new_nnode = 4;
            end
            
            valid_connec = obj.findConnecIndexes(coord_indexes,new_nnode);
            face_connec_raw = obj.mesh.connec(valid_connec,:);
            face_connec = obj.removeExtraNodes(face_connec_raw,coord_indexes,new_nnode);
        end
        
        function indexes = findConnecIndexes(obj,coord_indexes,nnode)
            number_of_valid_nodes_per_element = sum(ismember(obj.mesh.connec,coord_indexes),2);
            indexes = number_of_valid_nodes_per_element == nnode;
        end
    end
    
    methods (Static, Access = public)
        function quadrature = getQuadrature_Unfitted
            quadrature = Quadrature_Triangle;
        end
    end
    
    methods (Static, Access = private)
        function face_connec = removeExtraNodes(face_connec_raw,coord_indexes,nnode)
            valid_nodes = ismember(face_connec_raw,coord_indexes);
            
            face_connec = zeros(size(face_connec_raw,1),nnode);
            for i = 1:size(face_connec,1)
                face_connec(i,:) = face_connec_raw(i,valid_nodes(i,:));
            end
        end
    end
end

