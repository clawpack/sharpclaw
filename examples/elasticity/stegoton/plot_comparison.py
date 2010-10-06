from pyclaw.data import ClawPlotData
from pyclaw.plotting import plotframe

plotdata = ClawPlotData()
plotdata.outdir = '.'

# Figure:
plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)
plotfigure.kwargs = {'figsize':[5,3]}

# Axes:
plotaxes = plotfigure.new_plotaxes(name='Strain')
#plotaxes.xlim = [73,79]

plotitem = plotaxes.new_plotitem(name='SharpClaw 3600', plot_type='1d')
plotitem.plot_var = 0       # q[2] is the stress
plotitem.plotstyle = 's'
plotitem.color = 'b'   # could use 'r' or 'red' or '[1,0,0]'
plotitem.kwargs = {'linewidth':3,'markersize':10}

plotitem = plotaxes.new_plotitem(name='ClawPack 3600', plot_type='1d')
plotitem.outdir = '/users/ketch/research/claw42/fwave2/3600'
plotitem.plot_var = 0       # q[2] is the stress
plotitem.plotstyle = 'o'
plotitem.color = 'r'
plotitem.kwargs = {'linewidth':3,'markersize':10}

#plotitem = plotaxes.new_plotitem(name='ClawPack 28800', plot_type='1d')
#plotitem.outdir = '/users/ketch/research/claw42/fwave2/'
#plotitem.plot_var = 0       # q[2] is the stress
#plotitem.plotstyle = '-'
#plotitem.color = 'k'
#plotitem.kwargs = {'linewidth':3}


plotdata.plotframe(100)
