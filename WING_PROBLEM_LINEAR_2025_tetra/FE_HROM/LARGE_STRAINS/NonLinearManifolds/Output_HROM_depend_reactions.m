function  REACTIONS_time_HROM = Output_HROM_depend_reactions(DATA,BstRED_reactions,WeigthsECM,PK1STRESS)


nF = DATA.MESH.nstrain_F;
        REACTIONS_time_HROM = zeros(size(BstRED_reactions,2),size(PK1STRESS,2)) ;
        for itime = 1:size(WeigthsECM,2)
            w= WeigthsECM(:,itime) ;
            PoneST_w = zeros(size(PK1STRESS,1),1) ;
            for icomp = 1:nF
                icol = icomp:nF:size(PK1STRESS,1) ;
                PoneST_w(icol) = PK1STRESS(icol,itime).*w;
            end
            REACTIONS_time_HROM(:,itime) = BstRED_reactions'*PoneST_w ;
        end