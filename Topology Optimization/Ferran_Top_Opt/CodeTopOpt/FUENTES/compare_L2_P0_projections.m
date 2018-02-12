function [err1,err2] = compare_L2_P0_projections(phi,P_operator,Msmooth,A_nodal_2_gauss,coordinates,element,problembsc,dim)


F1 = faireF2(coordinates',element.conectivities',phi);

gamma1 = P_operator*F1; %more diffusive (big error), positive kernel, max principle holds
gamma2 = A_nodal_2_gauss*(Msmooth\F1); %middle diffusive (less error) , non positive kernel, max principle no holds
gamma_aver = cal_vol_mat(phi,dim,element,problembsc,coordinates)'; %less diffusion (no error), max principle holds

err1 = norm(gamma_aver - gamma1)/norm(gamma_aver);
err2 = norm(gamma_aver - gamma2)/norm(gamma_aver);