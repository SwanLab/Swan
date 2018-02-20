<<<<<<< HEAD
function [fixnodes,pnods] = compute_fix_nodes(BC,coordinates,lnodes,ftype,file_name)
pnods = [];fixnodes = [];

switch ftype
    case {'ELASTIC'}
        switch BC
            case 'PERIODIC'
                switch file_name
                    case {'RVE_Square','RVE_Square_Fine','RVE_Square_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_square(coordinates);
                    case {'RVE_Hexagonal','RVE_Hexagonal_Fine','RVE_Hexagonal_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_hexagonal(coordinates);
                    otherwise
                        error('RVE not detected');
                end
                
                % prescription of the corners
                ifix=0;
                for i=1:size(corners,2)
                    ifix=ifix+1;
                    fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=1; % idim
                    fixnodes(ifix,3)=0; % U_imp
                    ifix=ifix+1;
                    fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=2; % idim
                    fixnodes(ifix,3)=0; % U_imp
                end
                
            case 'NORMAL_DIRICLET'
                ifix=0;
                nfix = size(lnodes,1);
                for i=1:nfix
                    ifix=i;
                    fixnodes(ifix,1)=lnodes(i,1); % node
                    fixnodes(ifix,2)=lnodes(i,2); % idim
                    fixnodes(ifix,3)=0; % U_imp
                    %             ifix=ifix+1;
                    %             fixnodes(ifix,1)=lnodes(i); % node
                    %             fixnodes(ifix,2)=2; % idim
                    %             fixnodes(ifix,3)=0;
                end
                %         for i=1:size(coordinates,1);
                %             if (coordinates(i,2)<0.0025/2)
                %                 ifix=ifix+1;
                %                 fixnodes(ifix,1)=i;
                %                 fixnodes(ifix,2)=1; % idim
                %                 fixnodes(ifix,3)=0; % U_imp
                %             end
                %         end
                
                
             case 'NULL_DIRICLET'
                ifix=0;
                nfix = size(lnodes,1);
                for i=1:nfix
                    ifix=i;
                    fixnodes(ifix,1)=lnodes(i,1); % node
                    fixnodes(ifix,2)=lnodes(i,2); % idim
                    fixnodes(ifix,3)=lnodes(i,3); % U_imp
                    %             ifix=ifix+1;
                    %             fixnodes(ifix,1)=lnodes(i); % node
                    %             fixnodes(ifix,2)=2; % idim
                    %             fixnodes(ifix,3)=0;
                end
                %         for i=1:size(coordinates,1);
                %             if (coordinates(i,2)<0.0025/2)
                %                 ifix=ifix+1;
                %                 fixnodes(ifix,1)=i;
                %                 fixnodes(ifix,2)=1; % idim
                %                 fixnodes(ifix,3)=0; % U_imp
                %             end
                %         end   
                
            case 'SIMPLE_APOYO'
                ifix=1; i = 3278;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1; i = 4225;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1; i = 4225;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            case 'TODO_IMPUESTO'
                ifix=0;fixnodes=[];
                h=0.015625;
                %vertical izq
                for i=1:size(coordinates,1)
                    if (coordinates(i,1)<h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %vertical der
                for i=1:size(coordinates,1)
                    if (coordinates(i,1)>1-h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %horizontal inferior
                for i=1:size(coordinates,1)
                    if (coordinates(i,2)<h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %horizontal superior
                for i=1:size(coordinates,1)
                    if (coordinates(i,2)>1-h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
        end
        
    case {'THERMAL'}
        switch BC
            case 'NORMAL_DIRICLET'
                ifix=1;
                nfix = size(lnodes,1);
                for i=1:nfix
                    if lnodes(i,2) == 1
                        fixnodes(ifix,1)=lnodes(i,1); % node
                        fixnodes(ifix,2)=lnodes(i,2); % idim
                        fixnodes(ifix,3)=lnodes(i,3); % U_imp
                        ifix=ifix+1;
                    end
                end
                
            case 'NULL_DIRICLET'
                ifix=1;
                nfix = size(lnodes,1);
                for i=1:nfix
                    if lnodes(i,2) == 1
                        fixnodes(ifix,1)=lnodes(i,1); % node
                        fixnodes(ifix,2)=lnodes(i,2); % idim
                        fixnodes(ifix,3)= 0; % U_imp
                        ifix=ifix+1;
                    end
                end
                
                
            case 'PERIODIC'
                %fname = 'RVE04N3_PERIODIC';
                switch file_name
                    case {'RVE_Square','RVE_Square_Fine','RVE_Square_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_square(coordinates);
                    case {'RVE_Hexagonal','RVE_Hexagonal_Fine','RVE_Hexagonal_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_hexagonal(coordinates);
                    otherwise
                        error('RVE not detected');
                end
                
                ifix=1;
                for i=1:size(corners,2)
                   fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=1; % idim
                    fixnodes(ifix,3)=0; % U_imp
                    ifix=ifix+1;
                end
                
        end
end
end
=======
function [fixnodes,pnods] = compute_fix_nodes(BC,coordinates,dirichlet_data,ftype,file_name)
pnods = [];fixnodes = [];

switch ftype
    case {'ELASTIC'}
        switch BC
            case 'PERIODIC'
                switch file_name
                    case {'RVE_Square','RVE_Square_Fine','RVE_Square_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_square(coordinates);
                    case {'RVE_Hexagonal','RVE_Hexagonal_Fine','RVE_Hexagonal_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_hexagonal(coordinates);
                    otherwise
                        error('RVE not detected');
                end
                
                % prescription of the corners
                ifix=0;
                for i=1:size(corners,2)
                    ifix=ifix+1;
                    fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=1; % idim
                    fixnodes(ifix,3)=0; % U_imp
                    ifix=ifix+1;
                    fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=2; % idim
                    fixnodes(ifix,3)=0; % U_imp
                end
                
            case 'NORMAL_DIRICLET'
                ifix=0;
                nfix = size(dirichlet_data,1);
                for i=1:nfix
                    ifix=i;
                    fixnodes(ifix,1)=dirichlet_data(i,1); % node
                    fixnodes(ifix,2)=dirichlet_data(i,2); % idim
                    fixnodes(ifix,3)=0; % U_imp
                    %             ifix=ifix+1;
                    %             fixnodes(ifix,1)=dirichlet_data(i); % node
                    %             fixnodes(ifix,2)=2; % idim
                    %             fixnodes(ifix,3)=0;
                end
                %         for i=1:size(coordinates,1);
                %             if (coordinates(i,2)<0.0025/2)
                %                 ifix=ifix+1;
                %                 fixnodes(ifix,1)=i;
                %                 fixnodes(ifix,2)=1; % idim
                %                 fixnodes(ifix,3)=0; % U_imp
                %             end
                %         end
                
                
             case 'NULL_DIRICLET'
                ifix=0;
                nfix = size(dirichlet_data,1);
                for i=1:nfix
                    ifix=i;
                    fixnodes(ifix,1)=dirichlet_data(i,1); % node
                    fixnodes(ifix,2)=dirichlet_data(i,2); % idim
                    fixnodes(ifix,3)=dirichlet_data(i,3); % U_imp
                    %             ifix=ifix+1;
                    %             fixnodes(ifix,1)=dirichlet_data(i); % node
                    %             fixnodes(ifix,2)=2; % idim
                    %             fixnodes(ifix,3)=0;
                end
                %         for i=1:size(coordinates,1);
                %             if (coordinates(i,2)<0.0025/2)
                %                 ifix=ifix+1;
                %                 fixnodes(ifix,1)=i;
                %                 fixnodes(ifix,2)=1; % idim
                %                 fixnodes(ifix,3)=0; % U_imp
                %             end
                %         end   
                
            case 'SIMPLE_APOYO'
                ifix=1; i = 3278;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1; i = 4225;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1; i = 4225;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            case 'TODO_IMPUESTO'
                ifix=0;fixnodes=[];
                h=0.015625;
                %vertical izq
                for i=1:size(coordinates,1)
                    if (coordinates(i,1)<h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %vertical der
                for i=1:size(coordinates,1)
                    if (coordinates(i,1)>1-h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %horizontal inferior
                for i=1:size(coordinates,1)
                    if (coordinates(i,2)<h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
                %horizontal superior
                for i=1:size(coordinates,1)
                    if (coordinates(i,2)>1-h/5)
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=1; % idim
                        fixnodes(ifix,3)=0; % U_imp
                        ifix=ifix+1;
                        fixnodes(ifix,1)=i; % node
                        fixnodes(ifix,2)=2; % idim
                        fixnodes(ifix,3)=0; % U_imp
                    end
                end
        end
        
    case {'THERMAL'}
        switch BC
            case 'NORMAL_DIRICLET'
                ifix=1;
                nfix = size(dirichlet_data,1);
                for i=1:nfix
                    if dirichlet_data(i,2) == 1
                        fixnodes(ifix,1)=dirichlet_data(i,1); % node
                        fixnodes(ifix,2)=dirichlet_data(i,2); % idim
                        fixnodes(ifix,3)=dirichlet_data(i,3); % U_imp
                        ifix=ifix+1;
                    end
                end
                
            case 'NULL_DIRICLET'
                ifix=1;
                nfix = size(dirichlet_data,1);
                for i=1:nfix
                    if dirichlet_data(i,2) == 1
                        fixnodes(ifix,1)=dirichlet_data(i,1); % node
                        fixnodes(ifix,2)=dirichlet_data(i,2); % idim
                        fixnodes(ifix,3)= 0; % U_imp
                        ifix=ifix+1;
                    end
                end
                
                
            case 'PERIODIC'
                %fname = 'RVE04N3_PERIODIC';
                switch file_name
                    case {'RVE_Square','RVE_Square_Fine','RVE_Square_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_square(coordinates);
                    case {'RVE_Hexagonal','RVE_Hexagonal_Fine','RVE_Hexagonal_FineFine'}
                        [corners,pnods] = compute_periodic_boundary_nodes_hexagonal(coordinates);
                    otherwise
                        error('RVE not detected');
                end
                
                ifix=1;
                for i=1:size(corners,2)
                   fixnodes(ifix,1)=corners(i); % node
                    fixnodes(ifix,2)=1; % idim
                    fixnodes(ifix,3)=0; % U_imp
                    ifix=ifix+1;
                end
                
        end
end
end
>>>>>>> refs/remotes/origin/master
