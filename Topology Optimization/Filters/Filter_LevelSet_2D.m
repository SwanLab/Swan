classdef Filter_LevelSet_2D < Filter_LevelSet
    methods
        function obj = Filter_LevelSet_2D
            obj.max_subcells = 6;
            obj.nnodes_subelem = 3;
            obj.ndim = 2;
        end
        
        function getQuadrature_Unfitted(obj)
            obj.quadrature_unfitted = Quadrature_Triangle;
        end
        
        function getInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Triangle_Linear(obj.unfitted_mesh);
        end
        
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function M2=computeRHS_facet(obj,x,F)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeGlobalConnectivities;
%             obj.unfitted_mesh.plot;
            %             obj.unfitted_mesh.computeDvoluCut;
            
            [interpolation_facet,quadrature_facet] = obj.createFacet;
            interp_element = Interpolation.create(obj.mesh,obj.quadrature_fitted.order);
            
            facet_posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,interpolation_facet,quadrature_facet);
            
            shape_all = zeros(obj.nelem,obj.nnode);
            
            for ifacet = 1:size(obj.unfitted_mesh.connec,1)
                icell = obj.unfitted_mesh.cell_containing_subcell(ifacet);
                inode = obj.mesh.connec(icell,:);
                facet_posgp = facet_posgp_iso(:,:,ifacet);
                interp_element.computeShapeDeriv(facet_posgp');
                facet_deriv(:,:) = interpolation_facet.deriv(:,:,:);
                
                % !! How mapping is done for 2D cases??? !!
                t = [0; norm(diff(obj.unfitted_mesh.coord(obj.unfitted_mesh.connec(ifacet,:),:)))];
                dt_dxi = (facet_deriv'*t)/interpolation_facet.dvolu;
                
                f = (interp_element.shape*quadrature_facet.weigp')'*F(inode)/interpolation_facet.dvolu;
                shape_all(icell,:) = shape_all(icell,:) + (interp_element.shape*(dt_dxi.*quadrature_facet.weigp')*f)';
            end
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        function [P,active_nodes]=findCutPoints_Iso(obj,x,cut_elem,interpolation)
            gamma_1 = permute(x(obj.mesh.connec(cut_elem,:)),[2 3 1]);
            gamma_2 = permute([x(obj.mesh.connec(cut_elem,2:end)),x(obj.mesh.connec(cut_elem,1))],[2 3 1]);
            P1 = repmat(interpolation.pos_nodes,[1 1 size(cut_elem)]);
            P2 = repmat([interpolation.pos_nodes(2:end,:);interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [P,active_nodes]=findCutPoints_Global(obj,x,cut_elem,interpolation)
            index1 = permute(obj.mesh.connec(cut_elem,:),[2 3 1]);
            index2 = [permute(obj.mesh.connec(cut_elem,2:end),[2 3 1]);...
                permute(obj.mesh.connec(cut_elem,1),[2 3 1])];
            gamma_1=x(index1);
            gamma_2=x(index2);
            coord1 = obj.mesh.coord(:,1); coord2 = obj.mesh.coord(:,2);
            P1=[coord1(index1) coord2(index1)];
            P2=[coord1(index2) coord2(index2)];
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [interp_facet,quadrature_facet] = createFacet(obj)
            quadrature_facet = Quadrature.set('LINE');
            interp_facet = Line_Linear;
            quadrature_facet.computeQuadrature(obj.quadrature_fitted.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
    end
    
    methods (Static)
        function A=computeDvoluCut(elcrd)
            x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); x3 = elcrd(:,3,1); y3 = elcrd(:,3,2);
            A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
        
        function facets_connectivities = findFacetsConnectivities(elem_cutPoints_iso,interpolation,x,inode_global)
            if size(elem_cutPoints_iso,1) == 2
                facets_connectivities = [1 2];
            elseif size(elem_cutPoints_iso,1) == 4
                del_coord = [interpolation.pos_nodes;elem_cutPoints_iso];
                DT=delaunayTriangulation(del_coord);
                del_connec=DT.ConnectivityList;
                
                node_positive_iso = find(x(inode_global)>0);
                
                for idel = 1:length(node_positive_iso)
                    [a, ~] = find(del_connec==node_positive_iso(idel));
                    facets_connectivities(idel,:) = del_connec(a(end),del_connec(a(end),:)~=node_positive_iso(idel))-interpolation.nnode;
                end
            else
                error('Case still not implemented.')
            end
        end
        
        function djacob = mapping(elem_cutPoints_global,facets_connectivities,facet_deriv,dvolu)
            t = [0; norm(diff(elem_cutPoints_global(facets_connectivities,:)))];
            djacob = (facet_deriv'*t)/dvolu;
        end
    end
end

