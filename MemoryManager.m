classdef MemoryManager < handle
    
    methods (Access = public, Abstract)
        
        link(obj)
        allocateMemory(obj)
        freeSpareMemory(obj)
        transferData(obj)
        
    end
    
end

