function [gfunc] =  cal_nodgfunc_objective_macro_elast(igaus,gamma,element,problembsc,dim,strain,stres)
                                                       

% [chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void);
% [P_mas,P_menos] = polarization(igaus,element,problembsc,dim,chi,problembsc.flag_change_micro);

[P_interp] = derivative_constitutive_tensor(element,dim,gamma);

g = zeros(dim.nelem,1);

for istre=1:dim.nstre
     for jstre = 1:dim.nstre
        g = g + (stres(istre,:).*squeeze(P_interp(istre,jstre,:))'.*strain(jstre,:))'; 
     end
end

gfunc(igaus,:) = g;
%mean(g_mas./g_menos);


%std(g_mas./g_menos)

% gfunc(igaus,:)= -g_menos;
% outside= (phigp > 0);
% if any(outside)
%     gfunc(igaus,outside)= g_mas(igaus,outside);
% end

end 





