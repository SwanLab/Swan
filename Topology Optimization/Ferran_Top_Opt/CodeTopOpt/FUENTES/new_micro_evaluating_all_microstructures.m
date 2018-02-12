function [micr,Ener_Ant,Ener_Act] = new_micro_evaluating_all_microstructures(igrup,Group,stress,Ch_inv,micros_old,posgp,element,ndime,nnode,nelem,coordinates,ptype,weigp,nsnap,ngaus)

    elem_grup_log = Group(:,2) == igrup;
    elem_grup = Group(elem_grup_log,1);
                
    Ch_anter =  Ch_inv(:,:,micros_old(igrup));          
                
    Ener_Ant = 0;
    Ener = zeros(nsnap,1);
                
    for ielem_grup = 1:length(elem_grup)
        sigm_grup_elem = stress(:,:,elem_grup(ielem_grup));
        
        for igaus = 1:ngaus
            sim_grup_elem_gauss = sigm_grup_elem(igaus,:)';
            [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinates,coordinates,ptype);
            dvolu = weigp(igaus)*djacb(elem_grup(ielem_grup));
                       
            for i=1:3
                for j=1:3
                    Ener(:,1) = Ener(:,1) + (squeeze(sim_grup_elem_gauss(i)*Ch_inv(i,j,:)*sim_grup_elem_gauss(j)*dvolu));
                end
                
            end
               
             Ener_Ant = Ener_Ant + sim_grup_elem_gauss'*(Ch_anter)*sim_grup_elem_gauss*dvolu;
            
        end
        
    end
    
  
    [~,micr] = min(Ener(:,1));
    %[~,micr] = min(abs(0.9*Ener(igrup,micros_old(igrup)) - Ener(igrup,:)));
    Ener_Act = Ener(micr,1);



end