function Markers = MarkerMatrix(ntrial,ShowMarker,GiveMark,varargin);
% Different Markers
if  nargin == 0
    ntrial = 3;
    ShowMarker = 1;
elseif nargin == 1
    ShowMarker = 1 ;
elseif nargin == 2
    GiveMark = 'none';
end

% LOCAL INPUTS
ALEATORIO = 'YES';
% END LOCAL INPUTS

%% EXTRACTING INPUTS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varginOR = varargin ;
FdnamesInputs = {'ALEATORIO'};
AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
for id = 1:length(AuxDATA);
    eval(AuxDATA{id});
end



VecMark = {'diamond','pentagram','square','o','*','v','.','x','+',...
    '>','<','pentagram','hexagram','.'};

MarkersNew = cell(1,length(VecMark));
ppp = 1:length(VecMark);
switch ALEATORIO
    case 'YES'
        ppp = randperm(length(VecMark));
end


for i = 1:length(VecMark)
    MarkersNew{i} =  VecMark{ppp(i)};
end
VecMark =MarkersNew ;

if ShowMarker == 0
    VecMark = {'none'};
elseif ShowMarker == -1
    VecMark = {GiveMark} ;
end
nmark = length(VecMark);
iloc = 1;
nloc = 0 ;

for i = 1:ntrial
    if iloc>nmark
        iloc = 1 ;
        nloc = nloc+nmark ;
    end
    if i == 1
        Markers = {VecMark{iloc}} ;
    else
        Markers = cat(2,Markers,{VecMark{iloc}});
    end

    iloc = iloc + 1;

end


