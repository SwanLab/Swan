import time 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

tclock = dict()
def tic(ref=0):
    global tclock
    tclock[ref] = time.time()

def toc(ref=0):
    global tclock
    return format(time.time()-tclock[ref],"0.2f")+"s"

try:
    import colored as col

    def colored(text, color=None, attr=None):
        if color:
            text = col.stylize(text, col.fg(color))
        if attr:
            text = col.stylize(text, col.attr(attr))
        return text
except:
    def colored(text, color=None, attr=None):
        return text


    
class Display:
    """Text printing function and save in the log file self.log_file.

    Parameter:
        message  : the text to print
        level    : (default 0). The importance of the message (0: maximum). The message 
                   is displayed if level<=self.debug
    """

    def __init__(self, debug=0):
        self.debug = debug

    def __call__(self, message, level=0, debug=None, color=None, attr=None, end='\n',
                 flag=None):
        """ Display function with tunable level of verbosity

        INPUTS
        ------

        message        :   text to be printed
        level          :   level of importance of the message; will be actually
                           printed if debug >= level
        debug          :   current verbosity level
        color, attr    :   formattings with the `colored` package
        end            :   if set to '', will remove the final line carriage return
        flag           :   an extra indicator equal to None, 'stdout' or 'stderr',
                           the last two indicating that the text
                           passed to display comes from the standard output or
                           or error of a shell command.
                           Useful if display is overrided.
        """
        if debug is None:   
            debug = self.debug
        if color or attr:
            message = colored(message, color, attr)
        if debug >= level:
            print(message, end=end, flush=True)
    
_display = Display()
    
def display(message, level=0, debug=None, color=None, attr=None, end="\n", flag=None):  
    if color is None:   
        if level >= 1:
            color = 100-level
    _display(message, level, debug, color, attr, end, flag)

def display_iteration(it, J, G, H, x, level, debug, **kwargs):
    display('\n', level+1, debug)
    message = f'{it}. J='+format(J, '.6g')+' '
    message += 'G=['+",".join([format(x, '.4g') for x in G[:5]])
    if len(G)>5:
        message += "... ] "
        message += '(||G||=' + format(np.linalg.norm(G, np.inf), '.4g') + ') '
    else:
        message += "] "
    message += 'H=['+",".join([format(x, '.4g') for x in H[:5]])
    if len(H)>5:
        message += "... ] "
        message += 'max(H)=' + format(max(H), '.4g')
    else:
        message +="]"
    display(message, level, debug, **kwargs)
    if len(str(x)) < 100:
        display(f'x={x}', level+5, debug, **kwargs)

def set_plot_settings():
    mp.rcParams['contour.negative_linestyle'] = 'solid'
    mp.rcParams['lines.linewidth'] = 2
    mp.rcParams['axes.labelsize'] = 16
    mp.rcParams['xtick.labelsize'] = 16
    mp.rcParams['ytick.labelsize'] = 16
    mp.rcParams['legend.fontsize'] = 15
    mp.rcParams['figure.autolayout'] = True

def drawMuls(results, path_length=False, maxit=None, **kwargs):
    """Draw Lagrange multipliers of the ``results`` array   
    returned by the nlspace_solve function
        
    :param results: the array of results returned by :func:``nlspace_solve``.
    :param path_length: if set to `True`, then the abscissa of the plot is the  
                        length of the optimization path rather than the number of   
                        x
    :param maxit:       if set to an integer value ``N``, then the plot is truncated 
                        at iteration ``N``. 
    :param **kwargs**:  extra parameters that are supported by  
                        `matplotlib.pyplot.plot <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html>`_.
                            
    If ``results`` refer to an optimization problem with more than 10 constraints,  
    then only the first 10 multipliers are plot.
    """
    set_plot_settings()
    muls = list(zip(*results['muls']))
    if maxit is None:
        maxit = len(results['s'])
    title = kwargs.pop('title','Muls')
    linewidth = kwargs.pop('linewidth',2)
    if path_length:
        abscissa = results['s'][:maxit]
    else:
        abscissa = list(range(maxit))
    for i, mul in enumerate(muls):
        plt.plot(abscissa[:-1], mul[:maxit], linewidth=linewidth,
                 label=r'$\mu_'+str(i)+r'$'+f' - {title}', **kwargs)
        if i>10:
            break
    if path_length:
        plt.xlabel('s', fontsize=16)
    plt.title(title)
    plt.legend()
    plt.gca().legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
        
def drawHistories(results, start=0): 
    """Draw histories     
    object of the results of an optimization.   
    This function draws a `matplotlib.pyplot.subplots <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html>`_keys 
    object with the objective function (``results['J']``), the equality constraints (``results['G']``)  
    and the inequality constraints (``results['H']``) on the same plot.

    :param results: the array of results returned by :func:``nlspace_solve``.
    :param start: if set to non zero, the plots show iterations starting from `start`.
    """
    set_plot_settings()
    p = len(results['G'][0])
    q = len(results['H'][0])
    nplots = 1+p+q
    nrow = int(np.floor(np.sqrt(nplots)))
    ncol = int(np.ceil(float(nplots)/nrow))
    fig, axarray = plt.subplots(nrow, ncol)
    if nrow == 1 and ncol > 1:
        axarray = axarray[None, :]
    titles = ['J', *['G'+str(i) for i in range(p)],
              *['H'+str(i) for i in range(q)]]
    data = dict()
    iterations = range(start, len(results['J']))
    data = [results['J'], *list(map(list, zip(*results['G']))), *list(
        map(list, zip(*results['H'])))]
    for n in range(nplots):
        i = n // ncol
        j = n % ncol
        if nrow == 1 and ncol == 1:
            axarray.plot(iterations,
                         data[n][start:],label=self.config_name)
            axarray.set_title(titles[n])
            if len(cases)>1:
                axarray.legend()
        else:
            axarray[i, j].plot(iterations,
                               data[n][start:])
            axarray[i, j].set_title(titles[n])
    return fig, axarray

