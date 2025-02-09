function [c] = sumStruct(a,b)

theFields = fields(a);
for i = 1:length(theFields)
    c.(theFields{i}) = a.(theFields{i}) + b.(theFields{i});
end


end

