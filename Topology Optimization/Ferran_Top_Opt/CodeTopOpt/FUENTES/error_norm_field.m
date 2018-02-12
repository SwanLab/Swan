function enorm = error_norm_field(gp,gp0,Msmooth)

enodal = gp0 - gp;
enorm = (enodal'*Msmooth*enodal)/(gp'*Msmooth*gp);

end