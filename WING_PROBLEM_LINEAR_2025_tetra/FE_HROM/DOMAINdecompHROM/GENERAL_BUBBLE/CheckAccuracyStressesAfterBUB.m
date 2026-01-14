function CheckAccuracyStressesAfterBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
            INFO_RVE,BasisUdeform_bubble,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)


DATAoffline =DefaultField(DATAoffline,'CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES',0) ;
if DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES == 1
    
    if DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0        
        BasisUdeform_bubble.PhiDEFbs = PhiDEF ;
        BasisUdeform_bubble.PhiDEFcomp = GammaBUB;   
        DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES = 1 ;
        [~ ,~ ,~,~]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
            INFO_RVE,BasisUdeform_bubble,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    else
        error('Option not available (4-May-2025)')
    end
    
end