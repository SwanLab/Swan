function ClusterIdealCheck(DATAoffline,istep,iCL,qNEW1_all)


 if ~isempty(DATAoffline.IdealIndexCluster)  && ~isempty(qNEW1_all)
        % This is just for testing purposes 
        if  DATAoffline.IdealIndexCluster==1
            indMINcl = DATAoffline.IdealIndexCluster(istep) ;
        else            
            indIDEAL = DATAoffline.IdealIndexCluster(istep) ;
            %      if indIDEAL ~= indMINcl
            disp(['Current cluster = ',num2str(iCL),' DIST centr =',num2str(DIST_cluster(iCL)),' q1 =',num2str(qNEW1_all(iCL))])  ;
            disp(['Next cluster = ',num2str(indMINcl),' DIST centr =',num2str(DIST_cluster(indMINcl)),' q1 =',num2str(qNEW1_all(indMINcl))])  ;
            disp(['Ideal cluster =',num2str(indIDEAL),' DIST centr =',num2str(DIST_cluster(indIDEAL)),' q1 =',num2str(qNEW1_all(indIDEAL))])  ;
            %     end
        end
    end