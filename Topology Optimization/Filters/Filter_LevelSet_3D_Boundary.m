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
            M2 = obj.computeRHS(ones(size(x)));
            S = sum(M2);
            
            filter2D = Filter_P1_LevelSet_2D_Interior;
            %             filter2D.loadProblem(obj.diffReacProb.problemID,'MACRO');
            
            for idime = 1:obj.mesh.ndim
                for iside = 1:2
                    [face_mesh,valid_nodes] = obj.createFaceMesh(idime,iside);
                    unfitted_mesh2D = Mesh_Unfitted_2D_Interior(face_mesh,obj.interpolation_unfitted);
                    unfitted_mesh2D.computeMesh(x(valid_nodes));
                    
                    if ~isempty(unfitted_mesh2D.connec)
                        filter2D.setupFromMesh(face_mesh,'MACRO');
                        filter2D.preProcess;
                        filter2D.unfitted_mesh = unfitted_mesh2D;
                        M2 = filter2D.computeRHS(ones(size(x(valid_nodes))));
                        S = S + sum(M2);
                    end
                end
            end
        end
    end
    
    methods ( Access = private)
        function [face_mesh, valid_nodes] = createFaceMesh(obj,idime,iside)
            [face_coord,valid_nodes] = obj.getFaceCoordinates(idime,iside);
            %             face_connec = obj.getFaceConnectivities(face_coord,idime,iside);
            face_connec = obj.computeDelaunay(face_coord);
            face_mesh = Mesh;
            face_mesh = face_mesh.create(face_coord,face_connec);
        end
        
        function [face_coord, valid_coord] = getFaceCoordinates(obj,idime,iside)
            if iside == 1
                L = min(obj.mesh.coord(:,idime));
            else
                L = max(obj.mesh.coord(:,idime));
            end
            
            valid_coord = obj.mesh.coord(:,idime) == L;
            face_coord = obj.mesh.coord(valid_coord,:);
            face_coord = obj.removeExtraDimension(face_coord,idime);
            %             face_coord = unique(face_coord,'row');
        end
        
        %         function face_connec = getFaceConnectivities(obj,face_coord,idime,iside)
        %             indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(face_coord,obj.mesh.coord);
        %             nnode = 3;
        %             valid_cells = sum(ismember(obj.mesh.connec,indexes_in_global_matrix),2)==nnode;
        %             face_connec = obj.mesh.connec(valid_cells,:);
        %             %!! NOT WORKING !!
        %         end
        
        function indexes = findConnecIndexes(obj,coord_indexes,nnode)
            number_of_valid_nodes_per_element = sum(ismember(obj.mesh.connec,coord_indexes),2);
            indexes = number_of_valid_nodes_per_element == nnode;
        end
    end
    
    methods (Static, Access = public)     
        function djacob = mapping(points,dvolu)
            v1 = diff(points([1 2],:));
            v2 = diff(points([1 3],:));
            A = 0.5*norm(cross(v1,v2));
            djacob = A/dvolu;
        end
        
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
        
        function face_coord = removeExtraDimension(face_coord,idime)
            dimen = [1 2 3];
            face_coord = face_coord(:,dimen(dimen~=idime));
        end
        
        function connectivities = computeDelaunay(coordinates)
            DT = delaunayTriangulation(coordinates);
            connectivities = DT.ConnectivityList;
        end
    end
end

