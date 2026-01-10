function NODESlnNEW = RemovePointsLinesHEXA(NODESln,NODESpnt)

%dbstop('4')
if nargin == 0
    load('tmp.mat')
end

NODESlnNEW  = cell(size(NODESln))  ; 
%LINES = cell2mat(LINES) ; 

for i = 1:length(NODESln)
    
    
  
    
    
    NODESlnNEW{i}  = setdiff(  NODESln{i},NODESpnt) ; 
    
end

% Lines parallel to z-axis 

 
end

% 
% 
% NODESlnNEW = cell(size(NODESln)) ; 
%  for idimA = 1:size(NODESln,1)
%     for jdimA = 1:2
%         NODA = NODESln{idimA,jdimA} ;  % Plane A
%         CommonNodes = [] ;
%         for idimB =1:size(NODESln,1)
%             if  idimA~=idimB
%                 for jdimB = 1:2
%                     NODB = NODESln{idimB,jdimB} ;  % Plane B
%                     CommonNodes = [CommonNodes; intersect(NODA,NODB)];
%                     
%                 end
%             end
%         end
%         NODESlnNEW{idimA,jdimA} = setdiff(NODA,CommonNodes) ;
%         
%         
%     end
% end
% 
% end



%
% for iplane = 1:length(LINESf)
%     iname = LINESf{iplane} ;
%     for jplane = 1:length(LINESf)
%         jname = LINESf{jplane} ;
%         ijname = [iname,'_',jname] ;
%         iijj = intersect(NODESln.(iname),NODESln.(jname));
%         if  jplane>iplane && ~isempty(iijj) && length(iijj) ~= length(NODESln.(iname))
%
%             NODESpnt = setfield(NODESpnt,ijname,iijj) ;
%
%         end
%     end
% end