function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,OTHER_output] = DetermineBasis_elast_plast_EUCLID(SNAPdisp,DOFl,ind_elastic,ind_plastic,...
    DATAoffline,DATA_interp,OTHER_output)


DATAlocSVD.ISRELATIVE = 1;
TOLlocSVD = 1e-10 ;
[PhiMaster_lin,SSS,VVV] = SVDT(SNAPdisp(DOFl,ind_elastic),TOLlocSVD,DATAlocSVD)  ;
[PhiMaster_lin_plot,SSS,VVV] = SVDT(SNAPdisp(:,ind_elastic),TOLlocSVD,DATAlocSVD)  ;

% Inelastic part (orthogonal to PhiMaster_lin)
SNAPdisp_plast_L = SNAPdisp(DOFl,ind_plastic) - PhiMaster_lin*(PhiMaster_lin'*SNAPdisp(DOFl,ind_plastic)) ;


%TOLlocSVD = DATAoffline.errorDISP ;
[Phi_non_All,SSS,VVV] = SVDT(SNAPdisp_plast_L)  ;
nnn = norm(SNAPdisp(DOFl,:),'fro') ;


% Checking that the error is below prescribed tolerance
ERROR_phi = 1e20 ;
%nmodes = DATAoffline.nmodes_PLASTIC;
nmodes =1;
ERROR_phi_before = 1e40 ; 
while ERROR_phi > DATAoffline.errorDISP && ERROR_phi_before > ERROR_phi  && nmodes <= size(Phi_non_All,2)
    ERROR_phi_before = ERROR_phi; 
    Phi_non = Phi_non_All(:,1:nmodes) ;
    
    PhiALL = [PhiMaster_lin,Phi_non];
    ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*(PhiALL'* SNAPdisp(DOFl,:)) ;
    
    
    ERROR_phi = norm(ERROR_phi_mat,'fro')/nnn ;
    disp(['ERROR svd = ',num2str(ERROR_phi),' for ','nmodes PLASTIC =',num2str(nmodes)])
    
    nmodes = nmodes + 1;
    
end

Phi_non =  Phi_non_All(:,1:(nmodes-1)) ;



% 
% DATAoffline.nmodes_PLASTIC = 18 ; 

% % CHECKING WHETHER THE FIRST MODE IS INJECTIVE WITH RESPECT TO THE INPUT
% % PARAMETER (IMPOSED DISPLACEMENT)
% 
% figure(454)
% hold on
% xlabel('Snapshot')
% ylabel('Right singular vector')
% 
% for iii = 1:size(VVV,2)
%     plot(VVV(:,iii),'DisplayName',['S',num2str(iii), '=',num2str(SSS(iii)/SSS(1))]); 
% end
%  




% This is just for plotting purposes.
 
SNAPdisp_plast_plot= SNAPdisp(:,ind_plastic) - PhiMaster_lin_plot*(PhiMaster_lin_plot'*SNAPdisp(:,ind_plastic)) ; 

DATAsvd=[];
[Phi_non_plot,SSS,VVV] = SVDT(SNAPdisp_plast_plot)  ; 
Phi_non_plot = Phi_non_plot(:,1:size(Phi_non,2)) ; 
 
OTHER_output.Phi_To_Plot = [PhiMaster_lin_plot,Phi_non_plot] ; 


 
disp('***********************************************************')
disp(['Number of elastic displacement modes =',num2str(size(PhiMaster_lin,2)), ' '])
disp('***********************************************************')
 
disp('***********************************************************')
disp(['Number of inelastic displacement modes =',num2str(size(Phi_non,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
 
PhiMaster_nonl = Phi_non(:,1) ; 
PhiSlave_nonl =  Phi_non(:,2:end) ; 

% "Decoder function" 
%  d_L = PhiMaster_lin*qMASTER_lin +  PhiMaster_nonl*qMASTER_nonl +
%  PhiSlave_nonl*qSLAVE_nonl
%  where 
%   qSLAVE_nonl = f(qMASTER_nonl)
% Thus 
%  d_L = PhiMaster_lin*qMASTER_lin +  PhiMaster_nonl*qMASTER_nonl + PhiSlave_nonl* f(qMASTER_nonl)

qMASTER_nonl = PhiMaster_nonl'*SNAPdisp_plast_L;
qSLAVE_nonl = PhiSlave_nonl'*SNAPdisp_plast_L ; 
