function NODESplNEW = RemoveLinesPlanes(NODESpl,NODESln)

if nargin == 0
    load('tmp.mat')
end
NODESplNEW = cell(size(NODESpl)) ; 
 for idimA = 1:size(NODESpl,1)
    for jdimA = 1:2
        NODA = NODESpl{idimA,jdimA} ;  % Plane A
        CommonNodes = [] ;
        for idimB =1:size(NODESpl,1)
            if  idimA~=idimB
                for jdimB = 1:2
                    NODB = NODESpl{idimB,jdimB} ;  % Plane B
                    CommonNodes = [CommonNodes; intersect(NODA,NODB)];
                    
                end
            end
        end
        NODESplNEW{idimA,jdimA} = setdiff(NODA,CommonNodes) ;
        
        
    end
end

end



%
% for iplane = 1:length(PLANESf)
%     iname = PLANESf{iplane} ;
%     for jplane = 1:length(PLANESf)
%         jname = PLANESf{jplane} ;
%         ijname = [iname,'_',jname] ;
%         iijj = intersect(NODESpl.(iname),NODESpl.(jname));
%         if  jplane>iplane && ~isempty(iijj) && length(iijj) ~= length(NODESpl.(iname))
%
%             NODESln = setfield(NODESln,ijname,iijj) ;
%
%         end
%     end
% end