
# Set up the plots to be done for each frame

# This file is imported by the plotting routines and then
# setplot is called to set plot parameters specifying what is to be plotted.


#--------------------------
def setplot(plotdata):
#--------------------------

    plotdata.clearfigures()

    plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits=[-20.0,40]
    plotaxes.ylimits=[0.0,1.0]
    plotaxes.title = 'Solution'

    axis_xlim = [-20.0,40]
    axis_ylim = [0.0, 1.0]

    plotpoints = True
    plotline = True

    if plotpoints:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plotstyle = 'o'
        plotitem.color = 'r'
        plotitem.axis_xlim = axis_xlim
        plotitem.axis_ylim = axis_ylim
        plotitem.plot_var_title = 'Computed Solution'

    if plotline:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plotstyle = '-'
        plotitem.color = 'b'
        plotitem.axis_xlim = axis_xlim
        plotitem.axis_ylim = axis_ylim

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?


    return plotdata
