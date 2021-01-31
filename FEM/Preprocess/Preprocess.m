classdef Preprocess<handle
    % Only reads
    
    %% !! SUSCEPTIBLE OF BEING RENAMED !!
    
    properties
    end
    
    methods(Static)
        function data = readFromGiD(filename)
            
            run(filename)
            data = struct;
            if exist('gidcoord','var')
                coord=gidcoord;
                connec=gidlnods;
            end
            data.xpoints=coord;
            if length(data.xpoints(1,:))==3
                data.xpoints(:,4)=0;
            end
            if any(connec(:,length(connec(1,:))))==0
                data.connectivities=connec(:,1:(length(connec(1,:))-1));
            else
                data.connectivities=connec;
            end
            data.geometry = strjoin(Data_prb(1));
            data.problem_dim = strjoin(Data_prb(3));
            data.problem_type = strjoin(Data_prb(5));
            data.scale = strjoin(Data_prb(6));
            
            if strcmpi(data.problem_type,'elastic')
                if exist('dirichlet_data','var')
                    data.dirichlet_data = dirichlet_data;
                elseif exist('lnodes','var')
                    data.dirichlet_data = lnodes;
                else
                end
                
                if exist('pointload_complete','var')                
                    data.pointload = pointload_complete;
                else
                    data.pointload = [];
                end
            end
        end
        
        function [fixnodes,forces,full_dirichlet_data,Master_slave] = getBC_mechanics(filename)
            run(filename)
            if exist('lnodes','var')
                dirichlet_data=lnodes;
            end
            fixnodes{1,1} = dirichlet_data;
            
            if exist('pointload_complete','var')
                forces = pointload_complete;
            else
                forces = [];
            end
                        
            if exist('External_border_nodes','var')
                full_dirichlet_data= External_border_nodes;
            else
                full_dirichlet_data = [];
            end
            
            if ~isempty(full_dirichlet_data)
                full_dirichlet_data(:,2)=ones(length(full_dirichlet_data(:,1)),1);
                full_dirichlet_data(:,3)=zeros(length(full_dirichlet_data(:,1)),1);
            end
            
            if ~exist('Master_slave','var')
                Master_slave = [];
            end
            
        end
        
        function [fixnodes,forces,full_dirichlet_data,Master_slave] = getBC_fluids(filename,mesh,geometry,interp)
            run(filename)
            nelem= mesh.nelem;
            full_dirichlet_data=External_border_nodes;
            if ~isempty(full_dirichlet_data)
                full_dirichlet_data(:,2)=ones(length(full_dirichlet_data(:,1)),1);
                full_dirichlet_data(:,3)=zeros(length(full_dirichlet_data(:,1)),1);
            end
            
            if ~exist('Master_slave','var')
                Master_slave = [];
            end
            
            s.mesh = mesh;
            s.interpolation = interp{1};
            c = ConnecCoordFromInterpAndMesh(s);            
            c.compute();
            xpoints = c.coord;
            
            nnode   = length(xpoints(:,1));
            
            if (~isempty(velocity))
                ind=1;
                
                for inode = 1: nnode
                    if xpoints(inode,1) ==0  || xpoints(inode,1) == 1 ...
                            || xpoints(inode,2) == 0 || xpoints(inode,2) == 1
                        fixnodes(ind,:)=[inode 1 0];
                        fixnodes(ind+1,:)=[inode 2 0];
                        ind = ind+2;
                        
                    end
                end
                fixnodes_u=fixnodes;
                
            end
            
            if (~isempty(pressure))
                % if strcmp(geometry(2).interpolation.order,geometry(.order) ==1
                fixnodes_p = pressure;
                %                 else
                %                 ind2=1;
                %
                %                 for i = 1:length(pressure(:,1))
                %
                % %                         [ind,~] = find(interpolation_geometry.T == geometry_fixnodes_p(i,1)); % este para los que tienen minimo un nodo en la pared
                %                    % siguiente linea para los que tienen una cara en la pared
                %                 [ind,~] = find(external_elements(:,2:3) == geometry_fixnodes_p(i,1)); % troba l'element on esta imposada la p a la geometria
                %                 fixnodes_p(ind2: ind2+length(ind)-1,:) = [external_elements(ind,1) ones(length(ind),1)*geometry_fixnodes_p(i,2:3)]; % se li aplica la condicio al primer dels elements
                %                     ind2 = ind2+length(ind);
                %                 end
                %                 fixnodes_p = unique(fixnodes_p,'rows');
                %                 end
                %                 for i = 1:length(obj.fixnodes_p(:,1))
                %                     obj.iD_p(i) = obj.fixnodes_p(i,1)*nunkn_p - nunkn_p + obj.fixnodes_p(i,2);
                %                 end
            end
            
            if (~isempty(Vol_force))
                %                 ind=1;
                
                %                 for inode = 1:length(interpolation_variable(1).xpoints(:,1))
                %                     pos_node= num2cell(interpolation_variable(1).xpoints(inode,1:2));
                %                     f = cell2mat(force(pos_node{:}));
                %                     F(ind:ind + length(f)-1,:) = [[inode;inode] [1;2] f];
                %                     ind=ind+length(f);
                %                 end
                geom = geometry(1);
                
                quadrature = Quadrature.set(mesh.type);
                quadrature.computeQuadrature(interp{1}.order);

                geom.computeGeometry(quadrature,interp{1})
                xV = quadrature.posgp;
                xGauss = mesh.computeXgauss(xV);
                for ielem = 1:nelem
                    ind=1;
                    for igaus = 1:quadrature.ngaus
                        xG = xGauss(:,igaus,ielem);
                        pos_node= num2cell(xG);
                        f = cell2mat(Vol_force(pos_node{:}));
                        F(:,igaus,ielem) = f;
                        ind=ind+length(f);
                    end
                end
                
                forces = F;
            else
                forces = [];
            end
            clear fixnodes;
            fixnodes{1} = fixnodes_u;
            fixnodes{2} = fixnodes_p;
        end
        
        function forces_adjoint=getBC_adjoint(filename)
            run(filename)
            forces_adjoint = pointload_adjoint;
        end
    end
end
