function [deformation_all,phi_all,theta_all,num_of_snapshots,phi_matrix,theta_matrix,number_matrix,theta_step] = vademecum_spehere_mesh(nlevel,pint,phyisical_type)

switch phyisical_type
    case 'ELASTIC'
        close all
        cstring = 'rgbcmky';
        iteration = 1;
        phi_all = zeros(3,1);
        theta_all = zeros(3,1);
        
        phi_matrix = zeros(2^(nlevel-1) + 1,2^(nlevel-1) + 1);
        theta_matrix = zeros(2^(nlevel-1) + 1,1);
        number_matrix = zeros(2^(nlevel-1) + 1,2^(nlevel-1) + 1);
        theta_step = pi/2/(2^(nlevel-1));
        
        for ilevel = 1:nlevel
            
            idescr = 2^(ilevel-1) + 1;
            b = linspace(0,pi/2,idescr);
            c = linspace(0,pi/2,idescr);
            
            position = linspace(1,2.^([nlevel]-1)+1,2.^([ilevel]-1)+1);
            
            for iparalel = 1:idescr
                a_paralelo = acos(cos(b(iparalel))*cos(c(iparalel)));
                a = linspace(0,a_paralelo,iparalel);
                
                
                for imeridiano = 1:iparalel
                    A = acos((cos(a(imeridiano)) - cos(b(iparalel))*cos(c(iparalel)))/(1e-15 + sin(b(iparalel))*sin(c(iparalel))));
                    phi = - pi/4 + A;
                    theta = b(iparalel);
                    deformation = ([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]);
                    
                    
                    if ilevel == 1
                        phi_all(iteration,1) = phi;
                        theta_all(iteration,1) = theta;
                        deformation_all(iteration,:) = deformation;
                        
                        phi_matrix(position(imeridiano),position(iparalel)) = phi;
                        theta_matrix(iparalel) = theta;
                        number_matrix(position(imeridiano),position(iparalel)) = iteration;
                        
                        
                        %            %figure(1)
                        %            hold on
                        %           % plot(phi_all(iteration,1)*180/pi,theta_all(iteration,1)*180/pi,['+',cstring(mod(ilevel,7)+1)])
                        %            hold off
                        
                        iteration = iteration + 1;
                    end
                    
                    
                    
                    if ~((abs(phi_all - phi) <= 1e-6) & (abs(theta_all - theta) <= 1e-6) ) & ilevel > 1
                        
                        phi_all(iteration,1) = phi;
                        theta_all(iteration,1) = theta;
                        deformation_all(iteration,:) = deformation;
                        
                        phi_matrix(position(imeridiano),position(iparalel)) = phi;
                        theta_matrix(iparalel) = theta;
                        number_matrix(position(imeridiano),position(iparalel)) = iteration;
                        
                        
                        if pint
                            figure(1)
                            hold on
                            plot(phi_all(iteration,1)*180/pi,theta_all(iteration,1)*180/pi,['+',cstring(mod(ilevel,7)+1)])
                            hold off
                            
                        end
                        
                        
                        iteration = iteration + 1;
                    end
                    %spy(phi_matrix)
                    %pause(0.001)
                    
                end
                
            end
            num_of_snapshots(ilevel) = length(phi_all);
            
            if pint
                figure(19)
                %plot(num_of_snapshots)
                %[X,Y] = meshgrid(deformation_all(:,1),deformation_all(:,2));
                %Z = griddata(deformation_all(:,1),deformation_all(:,2),deformation_all(:,3),X,Y);
                %mesh(X,Y,Z)
                plot3(deformation_all(:,1),deformation_all(:,2),deformation_all(:,3),['o',cstring(mod(ilevel,7)+1)])
            end
            
        end
        
        
    case 'THERMAL'
        
        phi_all = [];
        for ilevel = 1:nlevel
            idescr = 2^(ilevel) + 1;
            if ilevel == 1
                phi_all = linspace(-pi/4,pi/4,idescr);
            else
                phi_all = [phi_all setdiff(linspace(-pi/4,pi/4,idescr),phi_all)];
            end
            num_of_snapshots(ilevel) = length(phi_all);
        end
        deformation_all = [cos(phi_all) sin(phi_all)];
        theta_all = zeros(size(phi_all));
        phi_matrix = 0;
        theta_matrix = 0;
        number_matrix = 0;
        theta_step = 0;


        
end






 
end
   