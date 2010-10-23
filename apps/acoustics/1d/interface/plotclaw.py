

def plotclaw():
    """
    Basic plotting script.  Execute from unix by
      $ python plotclaw.py
    Additional plot parameters are set in setplot.py.
    """
    
    from pyclaw.data import ClawPlotData
    from pyclaw.plotting import printframes

    plotdata = ClawPlotData()
    plotdata.setplot = 'setplot.py'   # module containing setplot function
                                      # to initialize plotting parameters

    # Parameters used only when creating html and/or pdf hardcopy
    # via pyclaw.plotting.printframes:
    plotdata.plotdir = 'plots'   # directory where plots will appear
    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_framesperline = 2         # layout of plots
    plotdata.latex_pdf = False               # also run pdflatex?

    # create html files of plots:
    printframes(plotdata, verbose=False)

#----------------------------------------------------------

if __name__=='__main__':
    # if executed at command line prompt, simply call the function:
    plotclaw()
