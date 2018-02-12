function [element,micros,Ener] = call_newmicros(stress,element,dim,problembsc,coordinates)
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; 

micros = element.micros;
ngaus = element.ngaus;
stress_vademecum = element.alpha0;

phi_matrix_total = element.phi_matrix_total;
theta_matrix_total = element.theta_matrix_total;
number_matrix_total = element.number_matrix_total;
theta_step = element.theta_step;



nstre = 3;

norm_stress =  sqrt(squeeze(stress(:,1,:).^2) + squeeze(stress(:,2,:).^2) + squeeze(stress(:,3,:).^2));

stress_unitary = stress;

for istres = 1:nstre
    stress_unitary(:,istres,:) = squeeze(stress(:,istres,:))./norm_stress;
end

   
for igaus = 1:ngaus
    for ielem = 1:nelem
        alpha = squeeze(stress_unitary(igaus,:,ielem))';
        alpha = sign(alpha(3))*alpha;
        [phi,theta] = s1mapping(alpha);
        itheta = floor(theta/(theta_step-1e-14))+1;
        
        if itheta == length(theta_matrix_total)
            itheta = itheta -1;
        end
        
        if itheta == 1
            interval_phi = 1:2;
            
        else
            interval_phi =  1:4*(itheta-1);
        end
        
            diff_i_phi = phi-phi_matrix_total( interval_phi,itheta);
            
            
            negative = interval_phi(diff_i_phi < 0);
             if isempty(negative)
                i_phi_plus_1_i_theta = 1;
            else
                i_phi_plus_1_i_theta = negative(1);
            end

            
            positive = interval_phi(diff_i_phi > 0);
            if isempty(positive)
                i_phi_i_theta = length(phi_matrix_total( interval_phi,itheta+1));
            else
                i_phi_i_theta = positive(end);
            end
         

            %[~,i_phi_i_theta] = min(abs(diff_i_phi(diff_i_phi < 0)));
            %[~,i_phi_plus_1_i_theta] = min(diff_i_phi(diff_i_phi > 0));
            interval_phi =  1:4*(itheta);
            diff_i_phi = phi-phi_matrix_total( interval_phi,itheta+1);
            
            negative = interval_phi(diff_i_phi < 0);
            if isempty(negative)
                i_phi_plus_1_i_theta_plus_1 = 1;
            else
                i_phi_plus_1_i_theta_plus_1 = negative(1);
            end
            

            positive = interval_phi(diff_i_phi > 0);
            if isempty(positive)
                i_phi_i_theta_plus_1 = length(phi_matrix_total( interval_phi,itheta+1));
            else
                i_phi_i_theta_plus_1 = positive(end);
            end
        
        indexes = [i_phi_i_theta                itheta;
            i_phi_plus_1_i_theta         itheta;
            i_phi_i_theta_plus_1         itheta+1;
            i_phi_plus_1_i_theta_plus_1  itheta+1];
        
        
        
        phi_all4 = [phi_matrix_total(indexes(1,1),indexes(1,2));phi_matrix_total(indexes(2,1),indexes(2,2));...
            phi_matrix_total(indexes(3,1),indexes(3,2));phi_matrix_total(indexes(4,1),indexes(4,2))];
        theta_all4 = [theta_matrix_total(itheta);theta_matrix_total(itheta);...
            theta_matrix_total(itheta+1);theta_matrix_total(itheta+1)];
        def_all4 = phitheta2strain(phi_all4,theta_all4);
        
        def_diff = def_all4 - ones(4,1)*alpha';
        %[~,imin] = min(sqrt(def_diff(:,1).^2) + squeeze(def_diff(:,2).^2) + squeeze(def_diff(:,3).^2));
        [~,imin] = min(sqrt((def_diff(:,1).^2) + (def_diff(:,2).^2) + (def_diff(:,3).^2)));
        

        
        new_micro(igaus,ielem) = number_matrix_total(indexes(imin,1),indexes(imin,2));
        error(igaus,ielem) = norm(element.alpha0(number_matrix_total(indexes(imin,1),indexes(imin,2)),:)-alpha');
        %   [~,i_phi_i_theta_plus_1] = min(abs(diff_i_phi(diff_i_phi < 0)));
        % [~,i_phi_plus_1_i_theta_plus_1] = min(diff_i_phi(diff_i_phi > 0));
        
        
        
        %
        %                 [phi_matrix_total(indexes(1,1),indexes(1,2)),phi_matrix_total(indexes(2,1),indexes(2,2))]
        %                 [phi_matrix_total(indexes(3,1),indexes(3,2)),phi_matrix_total(indexes(4,1),indexes(4,2))]
        %
        %                 phi
        %
        %                 [theta_matrix_total(itheta),theta_matrix_total(itheta+1)]
        %                 theta
        
        %plot3(def_all4(:,1),def_all4(:,2),def_all4(:,3),'+'), hold on,plot3(alpha(1),alpha(2),alpha(3),'r+')
        
        element.phi_gaus(igaus,ielem) = phi;
        element.theta_gaus(igaus,ielem) = theta;
    end
    
end

Ch = element.Ch0;
Ch_inv = element.Ch_inv0;
micros = new_micro;
for igaus = 1:element.ngaus
    element.Ch(igaus,:,:) = Ch(micros(igaus,:),:)';
    element.Ch_inv(igaus,:,:) = Ch_inv(micros(igaus,:),:)';
end    
    

element.micros = micros;



Energy = zeros(element.ngaus,dim.nelem);


for igauss = 1:element.ngaus
    
    [posgp,weigp] = cal_posgp_weigp(element.type,dim.ndime,dim.nnode,element.ngaus);
    [~,djacb] = cal_cartd(igauss,posgp,element,dim.ndime,dim.nnode,dim.nelem,coordinates,coordinates,problembsc.problemtype);
    dvolu = weigp(igauss)*djacb;
    
    Ch_inv =  vectorCh_2_tensorCh(squeeze(element.Ch_inv(igauss,:,:))');
    
    for i=1:3
        for j=1:3
            Energy(igauss,:) = Energy(igauss,:) + (squeeze(stress(igauss,i,:).*Ch_inv(i,j,:).*stress(igauss,j,:)).*dvolu)';
        end
        
    end
    
end

Ener = sum(Energy(:));


end