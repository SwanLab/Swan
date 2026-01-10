function selected_columnsLOC = SelectedComlumnsFun(DATAIN,iproject,nDOMx,nDOMy)

indX = DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject}{1} ;
        indY = DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject}{2} ;
        
          % JAHO-3-April-2019.  (Example)
          % NUMBER_OF_RVES =[5,5];  %  Number of slices in each test
          % DATA.RVES_SEL = {[2,3,4],[2,3,4]} ; 
          % Error detected. Columns were not properly selected 
          % Before:         selected_columnsLOC = sub2ind([nDOMx,nDOMy],indX,indY) ;
          % After : 
          iDOMxREP = repmat(indX(:),1,length(indY)) ; 
          iDOMyREP = repmat(indY(:)',length(indX),1) ;
           selected_columnsLOC = sub2ind([nDOMx,nDOMy],iDOMxREP(:),iDOMyREP(:)) ;