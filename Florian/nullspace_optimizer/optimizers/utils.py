class OptimizationResults():
    """ A wrapper class for dealing with optimization results   
        produced by null space. 
            
        :param Niterations: if set to a number ``N``, will save     
                            only the last ``N`` iterates of results['x']    
                                
        :param size_max:     maximum size for saving the constraints for iterations     
                             older than the last ``N`` steps."""
        
    def __init__(self, keys_group1, keys_group2,    
                 results=None, Niterations = None,  
                 size_max = 10, start=None):   
        self._results = dict()
        self._keys_group1 = keys_group1 
        self._keys_group2 = keys_group2
        self._keys = self._keys_group1 + self._keys_group2

        if results:
            self._results.update(results)
        if results is None or start==0:
            self._results = {k:[] for k in self._keys}
        if isinstance(start,int) and start>0:   
            for key in self._keys_group1:   
                self._results[key] = self._results[key][:start+1]
            for key in self._keys_group2:   
                self._results[key] = self._results[key][:start]

        self.checkResults()

        self._N = Niterations   
        self._SIZE_MAX = size_max
        

    def checkResults(self):
        """ Check the output dictionary in order to restart
            from the last iterate"""
        n = len(self._results[self._keys_group1[0]])
        for key in self._keys_group1:
            if key not in self._results:
                raise Exception(
                    f"Error, key {key} is missing in the results dictionary.")
            if len(self._results[key]) != n:
                raise Exception(f"Error, key {key} should have length n={n} "   
                                f" while it has length {len(self._results[key])}")
        for key in self._keys_group2:
            if key not in self._results:
                raise Exception(
                    f"Error, key {key} is missing in the results dictionary.")
            if len(self._results[key]) == n:
                self._results[key] = self._results[key][:-1]
            if n >0 and len(self._results[key]) != n-1:
                raise Exception(
                    f"Error, key {key} should have length n-1={n-1} while it has "  
                    f"length {len(self._results[key])}")
            if n == 0 and len(self._results[key]) != 0:
                raise Exception(
                    f"Error, key {key} should have length 0 while it has length "    
                    f"{len(self._results[key])}.")

    def save(self, key, value):
        self._results[key].append(value)
        if self._N and len(self._results[key])>self._N:  
            try:
                self._results[key][-self._N-1] = self._results[key][-self._N-1][:self._SIZE_MAX]
            except TypeError:   
                pass
            except IndexError:  
                pass
        if self._N and key=='x' and len(self._results[key])>self._N:
            self._results['x'][-self._N-1] = None


    def __getitem__(self, key): 
        if not key in self._keys:    
            raise Exception(f"Error, unknown key {key}. \"key\" should be one of"   
                            f"{self._keys}")
        return self._results[key]   
        
    def implementation(self):   
        return self._results
        
    def initialize(self):
        """ 
        Returns initialization and update ``results`` dictionary for restarting.    

        :returns: initial values for design variable ``x``, iteration number ``it`` 
                  ``s``, ``normxiJ_save``, ``muls``
        """
        if self._results[self._keys_group1[0]]:
            starting_values = { key : self._results[key][-1] for key in self._keys_group1 }
            for key in self._keys_group1:
                self._results[key] = self._results[key][:-1]
        else:   
            starting_values = None
        return starting_values
        

def check_params(default_parameters, params):
    if params is None:
        params = dict()
    for key in params:  
        if not key in default_parameters:   
            raise Exception("Error, parameter \""+key+"\" is not known. Admissible"
                            " parameters are \n" + "\n".join(default_parameters.keys()))
    for key in default_parameters:  
        params[key] = params.get(key,default_parameters[key])
    return params
