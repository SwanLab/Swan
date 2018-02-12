function  [element,micros,Ener_act,Ener_ant] = call_newmicros_by_group(stress,Group,element,dim,coordinates,ptype)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod;
Ch_inv = vectorCh_2_tensorCh(element.Ch_inv0);

[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,dim.ndime,dim.nnode,element.ngaus);

ngrup = max(Group(:,2));
nsnap =  size(Ch_inv,3);

micros_old = element.micros;
%Ener_Ant = zeros(ngrup,1);
%Ener_grup = zeros(ngrup,1);

micros = zeros(ngrup,1);
Ener_ant = 0;
Ener_act = 0;

%compactstress = compact_quadratic_stress(stress);
 
for igrup = 1:ngrup
    [micr,Ener_Ant_grup,Ener_Act_grup] = new_micro_evaluating_all_microstructures(igrup,Group,stress,Ch_inv,micros_old,posgp,element,ndime,nnode,nelem,coordinates,ptype,weigp,nsnap,ngaus);
  
   %[micr,Ener_Ant_grup,Ener_Act_grup] = new_micro_most_R6_parallel(igrup,Group,def_vadem,stress,Ch_inv,micros_old,wc,compactstress,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype,weigp,ngaus);
   
    micros(igrup) = micr;
    elem_grup_log = Group(:,2) == igrup;
    for igaus = 1:ngaus
    element.Ch(igaus,:,elem_grup_log) = repmat(element.Ch0(micr,:)',1,sum(elem_grup_log));
    element.micros(igaus,elem_grup_log) = micr;
    end
    Ener_ant = Ener_ant + Ener_Ant_grup;
    Ener_act = Ener_act + Ener_Act_grup;

end




end