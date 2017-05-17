import matplotlib.pyplot as p
import numpy as n
from matplotlib.patches import Polygon

def myax(fig_or_axes,
        conversion=None,rightaxlabel=None,
        autoscale=True,
        nudge=None,
        AL=.2,HL=.045,HW=.5,OH=.3,TW=.0035,PLW=0):
    '''
    myax(fig_or_axes,
        conversion=None,rightaxlabel=None,
        nudge=None,
        AL=.2,HL=.045,HW=.5,OH=.3,TW=.0035,PLW=0)
    Converts a figure with defined labels to K style.
    
    - Only required arg is axes hande or p.gca() (a figure handle is deprecated).
    -First two optional args are unit conversion function and the label if making a right axis
    -autoscale(=True):  If True, then the arrows are automatically scaled with the axes.
        If 'preserve', then the arrows are adjusted to be the exact same size as those on the standard 5x4 axes.
        If a float or int, then the aspect ratio is maintained but size adjusted.
        Both options hinge on the established defaults of HL=.045,HW=.5, and ,TW=.0035
        along with the mysty default axes size of 5x4.  So if these are ever messed with this could have issues.
        If False, then the defaults, or the values you input, are used.
    -nudge allows you to move the Y axis labels. By default they are position just outside the
     ticklabel bbox and at 3/4 in axes coordinates.  With nudge you specify a shift away from
     this default position as a fraction of the axis label's bbox
        Must be given as a list/tup/array specifying the shift as fraction of bbox
        i.e., nudge=('left',0.5,0.5) or nudge=(('left',0.5,0.5),('right',-0.5,0.5))
    -The rest are properties for the axis arrows.  These numbers are in Axes Coordinates.
    AL(=.2) is arrow length, including the head
    HL(=.045) is head length.
    HW(=.5) is head width as fraction of the head length!
    OH(=0.3) is overhang as fraction of head length!
    TW(=.0035) is tail width (in absolute axes values)
    PLW(=0) is the patch **kwarg linewidth
    Args for arrow patch pertain to the y-axis arrow, and are scaled acoordingly for the x-axis.
    The defaults are meant to optimize the appearance of the y-axis arrow on an axes with a 5x4 aspect ratio.
     '''
    from matplotlib.figure import Figure
    from matplotlib.axes._subplots import Subplot
    from matplotlib.axes._axes import Axes
    if type( fig_or_axes )  is Figure:
        fig = fig_or_axes
        ax = fig.gca()
    elif type( fig_or_axes ) in [Subplot,Axes]:
        ax = fig_or_axes
        fig = ax.get_figure()
    else:
        raise ValueError('Got invalid figure or axes type as first argument, {}'.format(fig_or_axes))
    
    # Test right now wether a second right axis is being called
    makeright = None not in [conversion,rightaxlabel]
    
    # Test right now for nudge
    nudgeleft, nudgeright = False, False
    if nudge is not None:
        nudge = n.asarray(nudge)
        if len(nudge.shape) == 1:
            if nudge[0] == 'left':
                nudge_left_x, nudge_left_y = nudge[1:].astype(float)
                nudgeleft = True
            elif nudge[0] == 'right':
                nudge_right_x, nudge_right_y = nudge[1:].astype(float)
                nudgeright = True
            else:
                raise ValueError('First arg of nudge must be "left" or "right"')
        elif len(nudge.shape) == 2:
            nudge_left_x, nudge_left_y = nudge[ nudge[:,0] == 'left',[1,2] ].astype(float).ravel()
            nudge_right_x, nudge_right_y = nudge[ nudge[:,0] == 'right',[1,2] ].astype(float).ravel()
            nudgeleft, nudgeright = True, True
        else:
            raise ValueError('Invalid Shape of args given to nudge')
    else:
         pass   
    
    # Figure dimensions and axis bounds in figure coordinates
    figwd,fight = fig.get_size_inches()
    bboxAX = fig.transFigure.inverted().transform(ax.bbox.get_points())
    x1,y1,x2,y2 = bboxAX.flatten()
    h2w = (y2-y1)*fight/(x2-x1)/figwd   #Ht to wt. conversion
    # Remove existing MyArrows; Currenty messed up if rightaxis already made
    for i in ax.get_children():
        if type(i) is MyArrow:
            i.remove()
    #An inverse transform that will take me from display to axes coordinates
    inv = ax.transAxes.inverted()
    #Handle for renderer
    R = fig.canvas.get_renderer()
    #Handle for x and y axis and labels
    axX = ax.xaxis
    xlab = axX.get_label()
    axY = ax.yaxis
    ylab = axY.get_label()
    # Turn off ticks on right and bottom
    ax.tick_params(right=False,top=False)
    
    # Check the labels to see whether they're mathtext fractions and increase size if so
    for lab in [xlab,ylab]:
        if lab.get_text().find(r'\frac') != -1:
            current_size = lab.get_size()
            lab.set_size( current_size * (30/18) )
            # 30 was found to be a good adjusted size for the incoming default of 18
    
    ########## X Label and Arrow ##########
    # Get the bounding box for the x tick labels and convert it to axes coords
    bboxXticklabs = inv.transform( axX.get_ticklabel_extents(R)[0].get_points() )
    xlab.set_va('top')
    xlab.set_horizontalalignment('left')
    # Move the axis label, and grab its new bbox
    axX.set_label_coords( (3/4), n.min(bboxXticklabs[:,1]), transform = ax.transAxes)
    bboxXlab = inv.transform(xlab.get_window_extent(R).get_points())
    # Set arrow parameters and draw the arrow
    left = n.min(bboxXlab[:,0])
    
    ### Autoscaling of arrow properties:
    if (autoscale in [True, 'preserve']) or (type(autoscale) in [int, float]):
        DY = (y2-y1)*fight
        DX = (x2-x1)*figwd
        if autoscale is 'preserve':
            HL=.045*(4/DY)
            HW=.5*5/DX*(DY/4)
            TW=.0035*5/DX
        else: 
            HL = .045*autoscale # n*True = n
            HW = .5*5/DX*(DY/4)
            TW = .0035*5/DX*(DY/4)*autoscale
    elif autoscale is False:
        pass
    else:
        raise ValueError('Invalid value "{}" given for argument autoscale.\nTrue, False, "preserve", or float/int required.'.format(autoscale))
     
    AL1 = AL            #(no *h2w because we want arrow length to scale with axis size)
    HL1=HL*h2w
    TW1=TW*(1/h2w)
    HW1=HW*(1/h2w)**2   # Squared because HW is a fraction of HL, and HL has been stretched by h2w
    OH1=OH              # No h2w here because OH is fraction of HL, which has already been scaled
    ar=MyArrow(left-AL1-HL1/2,n.mean(bboxXlab[:,1]),AL1,0,head_length=HL1,overhang=OH1,
                    head_width=HW1,tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
    ax.add_patch(ar)
    ar.set_clip_on(False)
    
    ########## Y label and Arrow #############
    # Get bbox for ytick labels
    bboxYticklabs = inv.transform( axY.get_ticklabel_extents(R)[0].get_points() )
    # Rotate and postion ylabel, then translate it based on the bbox overlap
    # We have to set its coords twice because rotating inevitibly leads to overlap with the ticklabels
    # So we set it once, calculate the overlap, then shift it once more to eliminate the overlap
    ylab.set_verticalalignment('center')
    ylab.set_horizontalalignment('center')
    ylab.set_rotation(0)
    xpos, ypos = n.min(bboxYticklabs[:,0]), (3/4) 
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes )
    bboxYlab = inv.transform(ylab.get_window_extent(R).get_points())
    xpos = xpos - ( n.max(bboxYlab[:,0]) - xpos )
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes)
    if nudgeleft:
        bbox_width, bbox_height = n.diff(bboxYlab, axis=0).ravel()
        xpos += bbox_width*nudge_left_x
        ypos += bbox_height*nudge_left_y
        axY.set_label_coords( xpos , ypos , transform = ax.transAxes)
        bboxYlab = inv.transform(ylab.get_window_extent(R).get_points())
    # Set arrow parameters and draw the arrow
    bot = n.min(bboxYlab[:,1])
    OH1=OH
    HL1=HL
    HW1=HW
    TW1=TW
    AL1=AL
    ar=MyArrow(xpos,bot-AL1-HL1/2,0,AL1,head_length=HL1,overhang=OH1, head_width=HW1,
                    tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
    ax.add_patch(ar)
    ar.set_clip_on(False)
    
    ###### Right Y-label #########
    ## Works the same way as previous y label formatting
    if makeright:
        axL, axR = ax, ax.twinx()
        # Doesn't quite work but is an improvement
        for i in axR.get_children():
            if type(i) is MyArrow:
                i.remove()
        y1, y2 = axL.get_ylim()
        axR.set_ylim(conversion(y1), conversion(y2))
        axR.set_ylabel(rightaxlabel, fontsize=ylab.get_size())
        axR.tick_params(labelsize=axY.get_ticklabels()[0].get_fontsize())
        axYR = axR.yaxis
        ylabR = axYR.get_label()
        bboxYRticklabs = inv.transform( axYR.get_ticklabel_extents(R)[1].get_points() )
        ylabR.set_verticalalignment('center')
        ylabR.set_horizontalalignment('center')
        ylabR.set_rotation(0)
        xpos, ypos = n.max(bboxYRticklabs[:,0]), (3/4)
        axYR.set_label_coords( xpos , ypos , transform = ax.transAxes )
        bboxYRlab = inv.transform(ylabR.get_window_extent(R).get_points())
        xpos = xpos + ( n.max(bboxYRlab[:,0]) - xpos )
        axYR.set_label_coords( xpos , ypos , transform = ax.transAxes)
        if nudgeright:
            bbox_width, bbox_height = n.diff(bboxYRlab, axis=0).ravel()
            xpos += bbox_width*nudge_right_x
            ypos += bbox_height*nudge_right_y
            axYR.set_label_coords( xpos , ypos ,     transform = ax.transAxes)
            bboxYRlab = inv.transform(ylabR.get_window_extent(R).get_points())
        bot = n.min(bboxYRlab[:,1])
        # Add arrow
        OH1=OH
        HL1=HL
        HW1=HW
        TW1=TW
        AL1=AL
        ar=MyArrow(xpos,bot-AL1-HL1/2,0,AL1,head_length=HL1,overhang=OH1, head_width=HW1,
                    tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
        axR.add_patch(ar)
        ar.set_clip_on(False)
        p.draw()
        return axL, axR
    else:
        p.draw()
        return ax

def makequad():
    ''' 
    Return a new figure with 2x2 tiled axes, properly spaced according to Fig.3
    in Chen, Scales, Kyriakides.
    Remember to call 'mysty-quad' style sheet so that fonts are scaled.
    Or, call 'mysty', and use rc_context to change font.size to 12 and axes.labelsize to 14.
    
    Returns: fig, ax00, ax01, ax10, ax11
    '''
    pad = 1
    hgap = 1.5
    vgap = 0.9
    axwt = 10/3
    axht = 8/3
    W = 2*pad + hgap + 2*axwt
    H = 2*pad + vgap + 2*axht
    fig1 = p.figure(figsize=(W,H))
    # Top left
    ax00 = fig1.add_axes([pad/W, (pad+vgap+axht)/H, axwt/W, axht/H])
    # Top right
    ax01 = fig1.add_axes([(pad+hgap+axwt)/W, (pad+vgap+axht)/H, axwt/W, axht/H])
    # Bottom left
    ax10 = fig1.add_axes([pad/W, pad/H, axwt/W, axht/H])
    # Bottom right
    ax11 = fig1.add_axes([(pad+hgap+axwt)/W, pad/H, axwt/W, axht/H])
    
    return fig1, ax00, ax01, ax10, ax11

def make21():
    '''
    Make a 2x1 subplot, properly layed out.
    Returns: fig, axtop, axbot
    Good with 'mysty.mplstyle'
    '''
    hpad = 1.5
    vpad = 1
    vgap = 1.5
    axht = 4
    axwt = 5
    W = 2*hpad+axwt
    H = 2*vpad+2*axht+vgap
    fig = p.figure(figsize=(W,H))
    axtop = fig.add_axes([hpad/W,(vpad+axht+vgap)/H,axwt/W,axht/H])
    axbot = fig.add_axes([hpad/W,vpad/H,axwt/W,axht/H])
    
    return fig, axtop, axbot
    
def make12():
    '''
    Make a 1x2 subplot, properly layed out.
    Returns: fig, axleft, axrt
    Good with 'mysty.mplstyle'
    '''
    hpad = 1.5
    vpad = 1
    hgap = 2
    axht = 4
    axwt = 5
    W = 2*hpad+2*axwt+hgap
    H = 2*vpad+axht
    fig = p.figure(figsize=(W,H))
    axleft = fig.add_axes([hpad/W,vpad/H,axwt/W,axht/H])
    axrt = fig.add_axes([(hpad+axwt+hgap)/W,vpad/H,axwt/W,axht/H])
    
    return fig, axleft, axrt    
    
def ksi2Mpa(ksi):
    """
    Multipls ksi by 6.89475908677537
    """
    return ksi * 6.89475908677537

def TickLabelFS(ax,fs=30,which='both'):
    '''
    Pass a figure handle.
    fs=fontsize
    Which specifies which axis (x or y) you want to change the ticklabel size of
    which = 'both' or two, 'x' or 0, 'y' or 1
    '''
    xls, yls = ax.xaxis.get_ticklabels(), ax.yaxis.get_ticklabels()
    if which in ['both' , 2]:
        itlist = [xls,yls]
    elif which in ['x',0]:
        itlist = [xls]
    elif which in ['y',1]:
        itlist = [yls]
    for labels in itlist:
        for i in labels:
            i.set_fontsize(fs)
    p.draw()
    return None

def data2axes(ax,coords):
    '''
    Converts data coords into axes coords.
    Most handy for adding annotation MyArrows within the axes box via ginput.
    '''
    return ax.transLimits.transform(coords)
    
def eztext(ax,text,loc='upper left',**kwargs):
    '''
    For placing text in the corners of the axis.
    Give the axis handle, the text, and the location.
    Location parameters:
        'upper left' = 'ul' = 0 = 'top left' = 'tl'
        'upper right' = 'ur' = 1 = 'top right' = 'tr'
        'lower right = 'lr'  = 2 = 'bottom right'  = 'br'
        'lower left' = 'll' = 3 = 'bottom left' = 'bl'
    '''
    if loc in ['upper left', 'ul', 0, 'top left', 'tl']:
        x, y = .01, .98
        ha, va = 'left', 'top'
    elif loc in ['upper right', 'ur', 1, 'top right', 'tr']:
        x, y = .99, .98
        ha, va = 'right', 'top'
    elif loc in ['lower right', 'lr', 2, 'bottom right', 'br']:
        x, y = .99, .01
        ha, va = 'right', 'bottom'
    elif loc in ['lower left', 'll', 3, 'bottom left', 'bl']:
        x, y = .01, .01
        ha, va = 'left', 'bottom'
    else:
        raise ValueError('Got unknown location, "{}"'.format(loc))
    tx = ax.text(x,y,text,ha=ha,va=va,transform=ax.transAxes,**kwargs)
    p.draw()
    return tx

def ezlegend(ax, loc='outside', nudge=(0,0), hl=0, txt2lc=True, **kwargs):
    ''' 
    A helper function to add a legend and format it K-style.
    Returns the legend handle.
    -First arg must be the axis handle.  Other args:
    -loc(='outside').
     If 'outside', then placed on the right, outside the axis frame.
     Otherwise, give it a MPL location (i.e., upper left).
    -nudge(=(0,0)), a tuple specifying in axes coords how much to nudge the location
     when 'outside' is the location.
    -hl(=0) is the line handlelength.  
     Set hl=None to use the default.
    -txt2lc(=True) controls whether or not to set text to the color of its corresponding line.
     If set to anything other than True, it will not do this.
    **kwargs:  Any other legend option can be passed.
    '''
    # Immediately check hl.  If it's nonzero, then we're leaving lines in and want to leave
    # the default handletextpad too.  Otherwise, handletextpad=0.
    if hl != 0:
        htp=None
    else:
        htp=0
    if loc == 'outside':
        leg = ax.legend(loc='center left',
                        bbox_to_anchor=(1.00 + nudge[0], .5 + nudge[1]),
                        bbox_transform=ax.transAxes,
                        handlelength=hl,
                        handletextpad=htp,  **kwargs)
    else:
        leg = ax.legend(loc=loc, handlelength=hl, handletextpad=htp, **kwargs)
    # Set texts to linecolor
    if txt2lc:
        [i.set_color(leg.get_lines()[k].get_color()) for k,i in enumerate(leg.get_texts())]
    # Set title size to same as texts
    p.setp(leg.get_title(), fontsize=leg.get_texts()[0].get_fontsize())
    # Return the legend
    return leg


def colorbar(ax,cbar,autoscale='preserve',AL=.2,HL=.045,HW=.5,OH=.3,TW=.0035,PLW=0):
    '''
    colorbar(ax,cbar,AL=.2,HL=.045,HW=.5,OH=.3,TW=.0035,PLW=0)
    
    *** Something seriously weird about how this (doesn't) works.
        When a script calling this function is just run, it renders 
        the arrow way far to the right of the cbar.
        When plotting interactively in ipython and a long code block is pasted into the QTconsole
        and f.colorbar() is contained in the code block, the arrow renders way far to the right.
        But if plotting interactively and I do everything the same except call f.colorbar, 
        then call f.colorbar on its own after the figure window has appeared, it works fine.
    ***
    
    Formats the colorbar label appropriately and adds an arrow that is identical to the 
    arrow on the figure axis.
    
    ax : the axis of the figure plot, NOT THE COLORBAR AXIS
    cbar : the colorbar, NOT THE COLORBAR AXIS
    The rest of the arguments are arrow props.
    Some day maybe I'll have it copy the y-axis arrow rather than specifying props.
    '''
    fig = ax.get_figure()
    cax = cbar.ax
    #An inverse transform that will take me from display to data coordinates
    # Notice that it is the plot axes, NOT the colorbar axes
    inv = ax.transAxes.inverted()
    #Handle for renderer
    R = fig.canvas.get_renderer()
    caxY = cax.yaxis
    ylab = caxY.get_label()
    
    if ylab.get_text().find(r'\frac') != -1:
            current_size = ylab.get_size()
            ylab.set_size( current_size * (30/18) )
    
    ticklabelbbox = inv.transform( caxY.get_ticklabel_extents(R)[1].get_points() )
    ylab.set_verticalalignment('center')
    ylab.set_horizontalalignment('center')
    ylab.set_rotation(0)
    xpos, ypos = n.max(ticklabelbbox[:,0]), (3/4)
    caxY.set_label_coords( xpos , ypos , transform = ax.transAxes )
    bboxLy = inv.transform(ylab.get_window_extent(R).get_points())
    xpos = xpos + ( n.max(bboxLy[:,0]) - xpos )
    caxY.set_label_coords( xpos , ypos , transform = ax.transAxes)
    # Set arrow parameters and draw the arrow
    bot = n.min(bboxLy[:,1])
  
    figwd,fight = fig.get_size_inches()
    bboxAX = fig.transFigure.inverted().transform(ax.bbox.get_points())
    x1,y1,x2,y2 = bboxAX.flatten()  
        ### Autoscaling of arrow properties:
    if (autoscale in [True, 'preserve']) or (type(autoscale) in [int, float]):
        DY = (y2-y1)*fight
        DX = (x2-x1)*figwd
        if autoscale is 'preserve':
            HL=.045*(4/DY)
            HW=.5*5/DX*(DY/4)
            TW=.0035*5/DX
        else: 
            HL = .045*autoscale # n*True = n
            HW = .5*5/DX*(DY/4)
            TW = .0035*5/DX*(DY/4)*autoscale
    elif autoscale is False:
        pass
    else:
        raise ValueError('Invalid value "{}" given for argument autoscale.\nTrue, False, "preserve", or float/int required.'.format(autoscale))
    
    
    OH1=OH
    HL1=HL
    HW1=HW
    TW1=TW
    AL1=AL
    ar=MyArrow(xpos,bot-AL1-HL1/2,0,AL1,head_length=HL1,overhang=OH1, head_width=HW1,
                    tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
    ax.add_patch(ar)
    ar.set_clip_on(False)
    p.draw()
    return ar

class MyArrow(Polygon):
    """
    Like Arrow, but lets you set head width and head height independently.
    """
    def __str__(self):
        return "MyArrow()"

    def __init__(self, x, y, dx, dy, tail_width=None, length_includes_head=True,
        head_width=None, head_length=None, shape='full', overhang=None,
        head_starts_at_zero=False, **kwargs):
        """
          *head_length*: Length of arrow head
            float or None. Default: 0.2 * total length
          *head_width*: Specified as fraction of head_length!
            float or None. Default: 0.5
           *overhang*: Specified as fraction of head legnth!
            Positive number means swept back, negative swept forward.
            float or None. Default: 0.2
          *tail_width*: width of full arrow tail
            float or None. Default: Length/50
          *shape*: ['full', 'left', 'right'] (default: 'full')
            draw the left-half, right-half, or full arrow
          *head_starts_at_zero*: [True | False] (default: False)
            if True, the head starts being drawn at coordinate 0
            instead of ending at coordinate 0.
          *length_includes_head*: [True | False] (default: False)
            True if head is to be counted in calculating the length.            
        Other valid kwargs (inherited from :class:`Patch`) are:
        %(Patch)s
        """
        distance = n.sqrt(dx ** 2 + dy ** 2)
        length = distance
        if head_length is None:
            head_length = .2 * length
        if not length_includes_head:
             length = distance + head_length
        if head_width is None:
            head_width = 0.6 
        if overhang is None:
            overhang = 0.3
        if tail_width is None:
            tail_width = .0035
        if not length:
            verts = []  # display nothing if empty
        else:
            hw, hl, oh, lw = head_width, head_length, overhang, tail_width
            left_half_arrow = n.array([
               [0.0, 0.0],                  # rear_center
                [0, -lw / 2.0],             # rear corner
                [length-hl, -lw / 2.0],     # feather-stem intersection
                [length-hl*(1+oh), -hw*hl / 2.0],          # sweep back
                [length, 0],                   # tip
                ])
            
            left_half_arrow += [-length, 0] # Keeps positioning consistent with MPL's arrow so I can use
                                            # their transformation below
            
            #p.figure(4)
            #p.plot(left_half_arrow[:,0],left_half_arrow[:,1],'-o')
            #for i in range(len(left_half_arrow)):
            #    p.text(left_half_arrow[i,0],left_half_arrow[i,1],'{}\n'.format(i),va='bottom',color='k',size=14)
                
            #if we're not including the head, shift up by head length
            if not length_includes_head:
                left_half_arrow += [head_length, 0]
            #if the head starts at 0, shift up by another head length
            if head_starts_at_zero:
                left_half_arrow += [head_length / 2.0, 0]
            #figure out the shape, and complete accordingly
            left_half_arrow = left_half_arrow[1:]
            if shape == 'left':
                coords = left_half_arrow
            else:
                right_half_arrow = n.flipud(left_half_arrow)[1:]*[1,-1]
                if shape == 'right':
                    coords = right_half_arrow
                elif shape == 'full':
                    # The half-arrows contain the midpoint of the stem,
                    # which we can omit from the full arrow. Including it
                    # twice caused a problem with xpdf.
                    coords = n.concatenate([left_half_arrow[:],
                                             right_half_arrow[:]])
                    coords = n.vstack( (coords,coords[0]) )[0:-1]
                    #p.figure(5)
                    #p.plot(coords[:,0],coords[:,1],'-o')
                    #for i in range(len(coords)):
                    #    #p.text(0,0,'what')
                    #    p.text(coords[i,0],coords[i,1],'{}\n'.format(i),va='bottom',color='k',size=14)
                else:
                    raise ValueError("Got unknown shape: %s" % shape)
            cx = float(dx) / distance
            sx = float(dy) / distance
            M = n.array([[cx, sx], [-sx, cx]])
            #verts = n.dot(coords, M)
            #p.figure(2)
            #p.plot(verts[:,0],verts[:,1],lw=1)
            #for i in verts:
            #    p.text(i[0],i[1],'({:.3f},{:.3f})'.format(i[0],i[1]),size=6)
            #verts += [x,y]
            #p.figure(6)
            #p.plot(verts[:,0],verts[:,1])
            #verts += [dx,dy]
            #p.figure(7)
            #p.plot(verts[:,0],verts[:,1])
            verts = n.dot(coords, M) + [x+dx, y+dy]
        Polygon.__init__(self, list(map(tuple, verts)), closed=True, **kwargs)    
