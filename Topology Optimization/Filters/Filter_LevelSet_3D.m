classdef Filter_LevelSet_3D < Filter_LevelSet
    properties % For Facets plotting
        surrounding_facets_connectivities
    end
    
    methods
        function obj = Filter_LevelSet_3D(unfitted_algorithm)
            obj@Filter_LevelSet(unfitted_algorithm);
            obj.max_subcells = 20;
            obj.nnodes_subelem = 4;
            obj.ndim = 3;
        end
        
        function getQuadratureDel(obj)
            obj.quadrature_del=Quadrature_Tetrahedra;
        end
        
        function mesh_del = getMeshDel(obj)
            mesh_del = obj.diffReacProb.mesh.duplicate;
            mesh_del.geometryType='TETRAHEDRA';
        end
        
        function getInterpolationDel(obj,mesh_del)
            obj.interp_del=Tetrahedra_Linear(mesh_del);
        end
        
        function setupUnfittedMesh(obj,x)
            switch obj.unfitted_mesh_algorithm
                case 'DELAUNAY'
                    obj.unfitted_mesh = Mesh_Unfitted_3D_Delaunay(obj.diffReacProb.mesh.duplicate,x,obj.diffReacProb.geometry.interpolation);
                case 'MARCHING_CUBES'
                    if strcmp(obj.geometry.type,'TETRAHEDRA')
                        obj.unfitted_mesh = Mesh_Unfitted_3D_MarchingCubes(obj.diffReacProb.mesh.duplicate,x,obj.diffReacProb.geometry.interpolation);
                    else
                        obj.unfitted_mesh = Mesh_Unfitted_3D_Delaunay(obj.diffReacProb.mesh.duplicate,x,obj.diffReacProb.geometry.interpolation);
                    end
            end
            obj.unfitted_mesh.findCutCells;
        end
        
        function M2 = computeRHS_facet(obj,x,F)
            [interp_facet,quadrature_facet] = obj.createFacet;
            interp_element = Interpolation.create(obj.diffReacProb.mesh,obj.quadrature.order);
            
            shape_all = zeros(obj.nelem,obj.nnode);
            [~,cut_elem]=obj.findCutElements(x,obj.connectivities);
            
            [P_iso,active_nodes_iso]=obj.findCutPoints_Iso(x,cut_elem,obj.geometry.interpolation);
            [P_global,active_nodes_global]=obj.findCutPoints_Global(x,cut_elem,obj.geometry.interpolation);
            
            %                         figure, hold on
            k = 0;
            for icut = 1:length(cut_elem)
                ielem = cut_elem(icut); inode_global = obj.connectivities(ielem,:);
                elem_cutPoints_iso = P_iso(active_nodes_iso(:,:,icut),:,icut);
                elem_cutPoints_global = P_global(active_nodes_global(:,:,icut),:,icut);
                
                interior_facets_local_connectivities = obj.computeFacetsLocalConnectivities(elem_cutPoints_iso,interp_element.pos_nodes,x(inode_global));
                
                for ifacet = 1:size(interior_facets_local_connectivities,1)
                    k = k+1;
                    facet_posgp = zeros(quadrature_facet.ngaus,interp_element.ndime);
                    %                     patch('vertices',facets_coordinates_global,'faces',facets_connectivities_global(k,:),'edgecolor',[0 1 0],...
                    %                         'facecolor','none','facelighting','phong')
                    %
                    for igaus = 1:quadrature_facet.ngaus
                        for idime = 1:interp_element.ndime
                            facet_posgp(igaus,idime) = interp_facet.shape(:,igaus)'*elem_cutPoints_iso(interior_facets_local_connectivities(ifacet,:),idime);
                        end
                    end
                    %                     plot3(facet_posgp(:,1),facet_posgp(:,2),facet_posgp(:,3),'xr')
                    
                    interp_element.computeShapeDeriv(facet_posgp');
                    facet_deriv(:,:) = interp_facet.deriv(:,:,:);
                    
                    djacob = obj.mapping(elem_cutPoints_global,interior_facets_local_connectivities(ifacet,:),facet_deriv,interp_facet.dvolu);
                    
                    f = (interp_element.shape*quadrature_facet.weigp')'*F(inode_global)/interp_facet.dvolu;
                    shape_all(ielem,:) = shape_all(ielem,:) + (interp_element.shape*(djacob.*quadrature_facet.weigp')*f)';
                    
                    %                     plot(obj.coordinates(obj.connectivities(ielem,:),1),obj.coordinates(obj.connectivities(ielem,:),2),'.-b'); plot(obj.coordinates(obj.connectivities(ielem,[1 4]),1),obj.coordinates(obj.connectivities(ielem,[1 4]),2),'.-b');
                    %                     plot(elem_cutPoints_global(facets_connectivities(i,:),1),elem_cutPoints_global(facets_connectivities(i,:),2),'-xr');
                    %                     title('Cut Elements & Cut points in GLOBAL coordinates'), axis('equal')
                end
            end
            %
            %             hold on
            %             patch('vertices',facets_coordinates_global,'faces',facets_connectivities_global,'edgecolor','none',...
            %                         'facecolor',[0 0 1],'facelighting','phong')
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        
        function [interior_facets_global_coordinates, interior_facets_global_connectivities] = computeInteriorFacets(obj,x)
            interp_element = Interpolation.create(obj.diffReacProb.mesh,obj.quadrature.order);
            
            [~,cut_elem]=obj.findCutElements(x,obj.connectivities);
            
            [P_iso,active_nodes_iso]=obj.findCutPoints_Iso(x,cut_elem,obj.geometry.interpolation);
            [P_global,active_nodes_global]=obj.findCutPoints_Global(x,cut_elem,obj.geometry.interpolation);
            
            interior_facets_global_coordinates = obj.findActiveCutPoints(P_global,active_nodes_global);
            
            interior_facets_global_connectivities = zeros(10*length(cut_elem),3);
            k = 0;
            for icut = 1:length(cut_elem)
                ielem = cut_elem(icut); inode_global = obj.connectivities(ielem,:);
                elem_cutPoints_iso = P_iso(active_nodes_iso(:,:,icut),:,icut);
                elem_cutPoints_global = P_global(active_nodes_global(:,:,icut),:,icut);
                
                interior_facets_local_connectivities = obj.computeFacetsLocalConnectivities(elem_cutPoints_iso,interp_element.pos_nodes,x(inode_global));
                
                for ifacet = 1:size(interior_facets_local_connectivities,1)
                    k = k+1;
                    interior_facets_global_connectivities(k,:) = obj.computeFromLocalToGlobalConnectivities(elem_cutPoints_global,interior_facets_global_coordinates,interior_facets_local_connectivities(ifacet,:));
                end
            end
            interior_facets_global_connectivities(k+1:end,:) = [];
        end
        
        function [P,active_nodes]=findCutPoints_Iso(obj,x,cut_elem,interpolation)
            iteration_1=interpolation.iteration(1,:);
            iteration_2=interpolation.iteration(2,:);
            gamma_1=permute(x(obj.connectivities(cut_elem,iteration_1)),[2 3 1]);
            gamma_2=permute(x(obj.connectivities(cut_elem,iteration_2)),[2 3 1]);
            P1=repmat(interpolation.pos_nodes(iteration_1,:),[1 1 size(cut_elem)]);
            P2=repmat(interpolation.pos_nodes(iteration_2,:),[1 1 size(cut_elem)]);
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes]=findCutPoints_Global(obj,x,cut_elem,interpolation)
            iteration_1=interpolation.iteration(1,:);
            iteration_2=interpolation.iteration(2,:);
            
            index1 = permute(obj.connectivities(cut_elem,iteration_1),[2 3 1]);
            index2 = permute(obj.connectivities(cut_elem,iteration_2),[2 3 1]);
            gamma_1=x(index1);
            gamma_2=x(index2);
            coord1 = obj.coordinates(:,1); coord2 = obj.coordinates(:,2); coord3 = obj.coordinates(:,3);
            P1=[coord1(index1) coord2(index1) coord3(index1)];
            P2=[coord1(index2) coord2(index2) coord3(index2)];
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [interp_facet,quadrature_facet] = createFacet(obj)
            quadrature_facet = Quadrature.set('TRIANGLE');
            interp_facet = Triangle_Linear(obj.diffReacProb.mesh);
            quadrature_facet.computeQuadrature(obj.quadrature.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
        
        function boundary_subfacets_connectivities = computeFacetsLocalConnectivities(obj,cutPoints_iso,pos_nodes, phi)
            del_coord = [pos_nodes; cutPoints_iso];
            DT=delaunayTriangulation(del_coord);
            subcells_connectivities=DT.ConnectivityList;
            
            interior_nodes = find(phi<=0); exterior_nodes = find(phi>0);
            
            % Find subcells formed by 1 interior node & 3 cutPoints
            counter_interior_nodes = obj.countInputNodesPerCell(subcells_connectivities,interior_nodes); %#ok<FNDSB>
            counter_exterior_nodes = obj.countInputNodesPerCell(subcells_connectivities,exterior_nodes); %#ok<FNDSB>
            boundary_subcells_connectivities = subcells_connectivities(counter_interior_nodes == 1 & counter_exterior_nodes == 0,:);
            
            boundary_subfacets_connectivities = zeros([size(boundary_subcells_connectivities,1),3]);
            for i = 1:size(boundary_subcells_connectivities,1)
                boundary_subfacets_connectivities(i,:) = boundary_subcells_connectivities(i,boundary_subcells_connectivities(i,:)>size(pos_nodes,1));
            end
            
            % !!!!!!!!!!!!!!!!!!!!!!! PLOTTING !!!!!!!!!!!!!!!!!!!!!!!!
            %             figure, hold on
            %             fac = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            %             patch('Faces',fac,'Vertices',[-1 -1 -1; -1 1 -1; 1 1 -1; 1 -1 -1; -1 -1 1; -1 1 1; 1 1 1; 1 -1 1],'FaceColor','w','FaceAlpha',0.0);
            %             for iconnec = 1:size(boundary_subfacets_connectivities,1)
            %                 plot3(del_coord(boundary_subfacets_connectivities(iconnec,:),1),del_coord(boundary_subfacets_connectivities(iconnec,:),2),del_coord(boundary_subfacets_connectivities(iconnec,:),3),'.-g')
            %                 plot3(del_coord(boundary_subfacets_connectivities(iconnec,[1 3]),1),del_coord(boundary_subfacets_connectivities(iconnec,[1 3]),2),del_coord(boundary_subfacets_connectivities(iconnec,[1 3]),3),'.-g')
            %             end
            %                                     patch('vertices',del_coord,'faces',boundary_subfacets_connectivities,'edgecolor',[0 0.5 0],...
            %                                                 'facecolor',[0 1 0],'facelighting','phong')
            
            %             plot3(cutPoints_iso(:,1),cutPoints_iso(:,2),cutPoints_iso(:,3),'og')
            %
            %             plot3(pos_nodes(phi>0,1),pos_nodes(phi>0,2),pos_nodes(phi>0,3),'+b')
            %             plot3(pos_nodes(phi<0,1),pos_nodes(phi<0,2),pos_nodes(phi<0,3),'<b')
            %
            %             axis equal
            %             view([115 20])
            %
            %             close
            boundary_subfacets_connectivities = boundary_subfacets_connectivities - size(pos_nodes,1);
        end
        
        function global_connectivities = computeFromLocalToGlobalConnectivities(obj,local_matrix_coordinates,global_matrix_coordinates,local_connectivities)
            indexes_in_global_matrix = obj.findCoordinatesIndexesInGlobalCoordinatesMatrix(local_matrix_coordinates,global_matrix_coordinates);
            global_connectivities = indexes_in_global_matrix(local_connectivities);
        end
        
        function [boundary_facets_coordinates,boundary_facets_connectivities] = computeBoundaryFacets(obj,x)
            [interior_facets_coordinates, interior_facets_connectivities] = obj.computeInteriorFacets(x);
            [surrounding_active_facets_coordinates,surrounding_active_facets_connectivities] = obj.computeSurroundingActiveFacets(x);
            
            boundary_facets_coordinates = [surrounding_active_facets_coordinates;interior_facets_coordinates];
            boundary_facets_connectivities = [surrounding_active_facets_connectivities;interior_facets_connectivities+size(surrounding_active_facets_coordinates,1)];
        end
        
        function computeSurroundingFacets(obj)
            surrounding_facets_coordinates_raw = zeros(size(obj.coordinates)); surrounding_facets_connectivities_raw = zeros(size(obj.connectivities,1),obj.ndim);
            k_coordinates = 0; k_connectivities = 0;
            for idime = 1:obj.ndim
                [surrounding_facets_coordinates_raw, surrounding_facets_connectivities_raw,k_coordinates,k_connectivities] = obj.computeBoxFaceCoordNConnec(surrounding_facets_coordinates_raw,surrounding_facets_connectivities_raw,idime,max(obj.coordinates(:,idime)),k_coordinates,k_connectivities);
                [surrounding_facets_coordinates_raw, surrounding_facets_connectivities_raw,k_coordinates,k_connectivities] = obj.computeBoxFaceCoordNConnec(surrounding_facets_coordinates_raw,surrounding_facets_connectivities_raw,idime,min(obj.coordinates(:,idime)),k_coordinates,k_connectivities);
            end
            surrounding_facets_coordinates_raw(k_coordinates+1:end,:) = [];
            surrounding_facets_connectivities_raw(k_connectivities+1:end,:) = [];
            
            obj.surrounding_facets_connectivities = obj.computeFromLocalToGlobalConnectivities(surrounding_facets_coordinates_raw,obj.coordinates,surrounding_facets_connectivities_raw);
            
            %             figure
            %             patch('vertices',obj.coordinates,'faces',obj.surrounding_facets_connectivities,...
            %                 'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
            %                 'facecolor','none','facelighting','flat')
            %             light
            %             axis equal off
        end
        
        function [surrounding_facets_coordinates, surrounding_facets_connectivities,k_coordinates,k_connectivities] = computeBoxFaceCoordNConnec(obj,surrounding_facets_coordinates,surrounding_facets_connectivities,idime,current_face_characteristic_coordinate,k_coordinates,k_connectivities)
            dimens = 1:obj.ndim;
            new_coordinates = obj.coordinates(obj.coordinates(:,idime) == current_face_characteristic_coordinate,:);
            DT = delaunayTriangulation(new_coordinates(:,dimens(dimens ~= idime)));
            new_connectivities = DT.ConnectivityList + k_coordinates;
            surrounding_facets_coordinates(k_coordinates+1:k_coordinates+size(new_coordinates,1),:) = new_coordinates;
            surrounding_facets_connectivities(k_connectivities+1:k_connectivities+size(new_connectivities,1),:) = new_connectivities;
            k_coordinates = k_coordinates + size(new_coordinates,1);
            k_connectivities = k_connectivities + size(new_connectivities,1);
        end
        
        function [surrounding_active_facets_coordinates,surrounding_active_facets_connectivities] = computeSurroundingActiveFacets(obj,x)
            [full_elem,cut_elem]=obj.findCutElements(x,obj.surrounding_facets_connectivities);
            surrounding_active_full_facets_connectivities = obj.surrounding_facets_connectivities(full_elem,:);
            surrounding_active_full_facets_coordinates = obj.coordinates;
            
            
            
            %             hold on

%             [interior_nodes_of_surrounding_active_cut_facets_coordinates,interior_nodes_of_surrounding_active_cut_facets_connectivities] = obj.computeSurfaceBoundaryInteriorSubcells(interp_element,obj.surrounding_facets_connectivities,x,cut_elem);
            
            %             plot3(interior_nodes_of_surrounding_active_cut_facets_coordinates(:,1),interior_nodes_of_surrounding_active_cut_facets_coordinates(:,2),interior_nodes_of_surrounding_active_cut_facets_coordinates(:,3),'.b')
            %
            %             surrounding_active_cut_facets_connectivities = obj.surrounding_facets_connectivities(cut_elem,:)
            
            interior_nodes_of_surrounding_active_cut_facets_coordinates = [];
            interior_nodes_of_surrounding_active_cut_facets_connectivities = [];
            
            surrounding_active_facets_coordinates = [surrounding_active_full_facets_coordinates;interior_nodes_of_surrounding_active_cut_facets_coordinates];
            surrounding_active_facets_connectivities = [surrounding_active_full_facets_connectivities;interior_nodes_of_surrounding_active_cut_facets_connectivities+size(surrounding_active_full_facets_coordinates,1)];
        end
    end
    
    methods (Static)
        function A=computeDvoluCut(elcrd)
            x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); z1=elcrd(:,1,3);
            x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); z2=elcrd(:,2,3);
            x3 = elcrd(:,3,1); y3 = elcrd(:,3,2); z3=elcrd(:,3,3);
            x4 = elcrd(:,4,1); y4 = elcrd(:,4,2); z4=elcrd(:,4,3);
            J=x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+...
                x4.*y1.*z2-x4.*y2.*z1-x1.*y3.*z4+x1.*y4.*z3+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1+x2.*y3.*z4-x2.*y4.*z3...
                -x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2;
            A=J/6;
        end
        
        function djacob = mapping(elem_cutPoints_global,facets_connectivities,facet_deriv,dvolu)
            % !! Remove facet_deriv --> no longer used !!
            v1 = diff(elem_cutPoints_global(facets_connectivities([1 2]),:));
            v2 = diff(elem_cutPoints_global(facets_connectivities([1 3]),:));
            A = 0.5*norm(cross(v1,v2));
            djacob = A/dvolu;
        end
        
        function all_cutPoints_global = findActiveCutPoints(P_global,active_nodes_global)
            P_global_x = P_global(:,1,:); P_global_y = P_global(:,2,:); P_global_z = P_global(:,3,:);
            all_cutPoints_global(:,1) = P_global_x(active_nodes_global); all_cutPoints_global(:,2) = P_global_y(active_nodes_global); all_cutPoints_global(:,3) = P_global_z(active_nodes_global);
            
            all_cutPoints_global = unique(all_cutPoints_global,'rows','stable');
        end
        
        function indexes_in_global_matrix = findCoordinatesIndexesInGlobalCoordinatesMatrix(coordinates_local,coordinates_global)
            indexes_in_global_matrix = zeros(1,size(coordinates_local,1));
            for inode = 1:size(coordinates_local,1)
                match = true(size(coordinates_global,1),1);
                for idime = 1:size(coordinates_local,2)
                    match = match & coordinates_global(:,idime) == coordinates_local(inode,idime);
                end
                indexes_in_global_matrix(inode) = find(match,1);
            end
        end
        
        function counter = countInputNodesPerCell(connectivities,nodes)
            counter = zeros(size(connectivities,1),1);
            for iconnec = 1:size(nodes,1)
                match = false(size(connectivities,1),1);
                for inode = 1:size(connectivities,2)
                    match = match | connectivities(:,inode) == nodes(iconnec);
                end
                counter = counter + match ;
            end
        end
    end
end

