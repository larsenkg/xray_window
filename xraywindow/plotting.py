import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd

ELEMENTS            = pd.DataFrame(columns=["Element", "Energy"])
ELEMENTS["Element"] = pd.Series(["Li", "Be", "B", "C", "N", "O", "Si"])
ELEMENTS["Energy"]  = pd.Series([54.33, 108.5, 183.3, 277, 392.4, 524.9, 1739.9])
# ELEMENTS["Element"] = pd.Series(["Li", "Be", "B", "C", "N", "O", "F", "Na", "Mg", "Al", "Si"])
# ELEMENTS["Energy"]  = pd.Series([54.33, 108.5, 183.3, 277, 392.4, 524.9, 676.8, 1040.98, 1253.6, 1486.7, 1739.9])
ELEMENTS            = ELEMENTS.set_index("Element")

# Add vertical lines for common elements
def add_elem_lines(
    ax         = None, 
    elements   = ELEMENTS, 
    colors     = [0,0,0,0.4], 
    ymin       = -1, 
    ymax       = 1, 
    min_energy = 0, 
    max_energy = 1200, 
    fontsize   = 8
):
    if ax is None:
        ax = plt.gca()
        
    query_str = f"Energy > {min_energy} and Energy < {max_energy}"
    
    ax.vlines("Energy", ymin, ymax, colors=colors, data=elements.query(query_str), label=None)
    
    for e in elements.query(query_str).iterrows():
        ax.text(e[1]["Energy"], ymax*1.03, e[0], rotation="vertical", horizontalalignment="center", fontsize=fontsize)
        
    return ax
        
def plot_transmission(
    x          = "Energy", 
    y          = "Transmission", 
    data       = None, 
    labels     = None,
    colors     = None,
    linestyles = None,
    add_lines  = True, 
    use_log    = False,
    use_ylog   = False,
    use_perc   = False,
    ymin       = 0, 
    ymax       = 1, 
    min_energy = 0, 
    max_energy = 1200,
    use_grid   = True,
    ax         = None,
    line_kws   = {}
):
    """
    Plot one or more transmission spectra
    
    Parameters
    ----------
    x, y : key in data
        Keys referencing columns in data.
    
    data : pd.DataFrame or list of pd.DataFrame objects
        Single or multiple pd.DataFrames. Will plot same columns from each.
        
    labels : list of str
        List of strings containing labels for data. Must be in the same order as data.
        
    
    """
    if type(data) == pd.DataFrame:
        data = [data]
        
    if ax is None:
        ax = plt.gca()
        
    if labels is None:
        labels = [""]*len(data)
        
    if colors is None:
        colors = [None]*len(data)
        
    if linestyles is None:
        linestyles = [None]*len(data)
                
    for d, label, color, linestyle in zip(data, labels, colors, linestyles):
        ax.plot(x, y, data=d, label=label, color=color, linestyle=linestyle, **line_kws)

    
    
    
    # Adjust ymin so 0% is not right at the bottom
    if ymin is None:
        ymin, _ = ax.get_ylim()
    if ymax is None:
        _, ymax = ax.get_ylim()
    
    if not use_ylog:
        ymin = ymin - ymax*0.05
    
    ax.legend()
    
    # Add vertical lines for common elements
    if add_lines:
        if use_ylog:
            ymax_tmp = 1.75*ymax
        else:
            ymax_tmp = ymax
        ax = add_elem_lines(ymin=ymin, ymax=ymax_tmp, min_energy=min_energy, max_energy=max_energy, ax=ax)

    ax.set_xlabel("Energy (eV)")
    if use_perc:
        ax.set_ylabel("Transmission (%)")
    else:
        ax.set_ylabel("Transmission")
        
    if use_log:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x,pos: f"{x:.0f}"))
    if use_ylog:
        ax.set_yscale('log')
    if use_perc:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(ymax))
        
    ax.set_xlim(min_energy, max_energy)
    ax.set_ylim(ymin, ymax)
    
    if use_grid:
        ax.grid(axis='y')
    
    plt.tight_layout()
    return ax