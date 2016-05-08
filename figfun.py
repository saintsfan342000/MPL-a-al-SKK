import matplotlib.pyplot as p
import numpy as n
from matplotlib.patches import Polygon

def myax(fig,conversion=None,rightaxlabel=None,AL=.22,HL=.045,OH=.3,TW=.0045,HW=.5,PLW=0):
    '''
    Converts a figure with to K style.
    Requires figure hande or p.gcf()
    Optional args are unit conversion on label if doing a second right ax,
    then properties for the axis arrows.  These numbers are in Axes Coordinates.
    AL(=.22) is arrow length, including the head
    HL(=.045) is head length.
    HW(=.5) is head width as fraction of the head length!
    OH(=2.9) is overhang as fraction of head length!
    TW(=.0045) is tail width (in absolute axes values)
    PLW(=0) is the patch **kwarg linewidth
    Args for arrow patch pertain to the y-axis arrow, and are scaled acoordingly for the x-axis.
    The defaults are meant to optimize the appearance of the y-axis arrow.
    If conversion and rightaxlabel are both not none, then the twinx() method is used to create a right ax
     '''
    #Handle for figure axes
    ax = fig.gca()
    #Test right now wether a second right axis is being called
    makeright = None not in [conversion,rightaxlabel]
    # Figure dimensions and axes bounds points in figure coordinates
    figwd,fight = fig.get_size_inches()
    bboxAX = fig.transFigure.inverted().transform(ax.bbox.get_points())
    x1,y1,x2,y2 = bboxAX.flatten()
    h2w = (y2-y1)*fight/(x2-x1)/figwd   #Ht to wt. conversion
    for i in ax.get_children():
        if type(i) is MyArrow:
            i.remove()
    #An inverse transform that will take me from display to data coordinates
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
            lab.set_size(30)

    ########## X Label and Arrow ##########
    # Get the bounding box for the x tick labels and convert it to axes coords, and set position accordingly
    bboxX = inv.transform( axX.get_ticklabel_extents(R)[0].get_points() )
    xlab.set_va('top')
    xlab.set_horizontalalignment('left')
    axX.set_label_coords( (3/4), n.min(bboxX[:,1]), transform = ax.transAxes)
    bboxLx = inv.transform(xlab.get_window_extent(R).get_points())
    left = n.min(bboxLx[:,0]) 
    AL1 = AL            #(no *h2w because we want arrow length to scale with axis size)
    HL1=HL*h2w
    TW1=TW*(1/h2w)
    HW1=HW*(1/h2w)**2   # Squared because HW is a fraction of HL, and HL has been stretched by h2w
    OH1=OH              # No h2w here because OH is fraction of HL, which has already been scaled
    ar=MyArrow(left-AL1-HL1/2,n.mean(bboxLx[:,1]),AL1,0,head_length=HL1,overhang=OH1,
                    head_width=HW1,tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
    ax.add_patch(ar)
    ar.set_clip_on(False)
    
    ########## Y label #############
    # Get bbox for ytick labels
    bboxY = inv.transform( axY.get_ticklabel_extents(R)[0].get_points() )
    # Rotate and postion ylabel, then translate it based on the bbox overlap
    ylab.set_verticalalignment('center')
    ylab.set_horizontalalignment('center')
    ylab.set_rotation(0)
    #xpos, ypos = n.min(bboxY[:,0]), (3/4) * n.abs(n.diff(bboxY[:,1])[0])
    xpos, ypos = n.min(bboxY[:,0]), (3/4) 
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes )
    bboxLy = inv.transform(ylab.get_window_extent(R).get_points())
    xpos = xpos - ( n.max(bboxLy[:,0]) - xpos )
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes)
    bot = n.min(bboxLy[:,1])
    # Add arrow
    OH1=OH
    HL1=HL
    HW1=HW
    TW1=TW
    AL1=AL
    ar=MyArrow(xpos,bot-AL1,0,AL1,head_length=HL1,overhang=OH1, head_width=HW1,
                    tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
    ax.add_patch(ar)
    ar.set_clip_on(False)
    
    ###### Right Y-label #########
    if makeright:
        axL, axR = ax, ax.twinx()
        y1, y2 = axL.get_ylim()
        axR.set_ylim(conversion(y1), conversion(y2))
        axR.set_ylabel(rightaxlabel)
        axYR = axR.yaxis
        ylabR=axYR.get_label()
        bboxYR = inv.transform( axYR.get_ticklabel_extents(R)[1].get_points() )
        ylabR.set_verticalalignment('center')
        ylabR.set_horizontalalignment('center')
        ylabR.set_rotation(0)
        xpos, ypos = n.max(bboxYR[:,0]), (3/4)
        axYR.set_label_coords( xpos , ypos , transform = ax.transAxes )
        bboxLyR = inv.transform(ylabR.get_window_extent(R).get_points())
        xpos = xpos + ( n.max(bboxLyR[:,0]) - xpos )
        axYR.set_label_coords( xpos , ypos , transform = ax.transAxes)
        bot = n.min(bboxLyR[:,1])
        # Add arrow
        OH1=OH
        HL1=HL
        HW1=HW
        TW1=TW
        AL1=AL
        ar=MyArrow(xpos,bot-AL1,0,AL1,head_length=HL1,overhang=OH1, head_width=HW1,
                    tail_width=TW1,lw=PLW,transform=ax.transAxes,color='k')
        axR.add_patch(ar)
        ar.set_clip_on(False)

    p.draw()
    return None

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
        
        #left-AL1,n.mean(bboxLx[:,1]),AL1,0,overhang=OH1,transform=ax.transAxes,
        #            color='k',head_length=HL1,length_includes_head=True,
        #            head_width=HW1,tail_width=TW1,lw=PLW,rasterized=rast
        
        
        
        distance = n.sqrt(dx ** 2 + dy ** 2)
        length = distance
        if head_length is None:
            print('head_length is None')
            head_length = .2 * length
        if not length_includes_head:
             length = distance + head_length
        if head_width is None:
            print('head_width is None')
            head_width = 0.6 
        if overhang is None:
            print('overhang is None')
            overhang = 0.3
        if tail_width is None:
            print('tail_width is None')
            tail_width = length / 50
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

def ksi2Mpa(ksi):
    """
    Multipls ksi by 6.89475908677537
    """
    return ksi * 6.89475908677537

def TickLabelFS(fig,fs=30,which='both'):
    '''
    Pass a figure handle.
    fs=fontsize
    Which specifies which axis (x or y) you want to change the ticklabel size of
    which = 'both' or two, 'x' or 0, 'y' or 1
    '''
    xls, yls = fig.gca().xaxis.get_ticklabels(), fig.gca().yaxis.get_ticklabels()
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
    return ax.transAxes.inverted().transform( ax.transData.transform( coords ) )

class FancyArrow_original(Polygon):
    """
    Like Arrow, but lets you set head width and head height independently.
    """

    def __str__(self):
        return "FancyArrow()"

    def __init__(self, x, y, dx, dy, width=0.001, length_includes_head=False,
        head_width=None, head_length=None, shape='full', overhang=0,
        head_starts_at_zero=False, **kwargs):
        """
        Constructor arguments
          *width*: float (default: 0.001)
            width of full arrow tail

          *length_includes_head*: [True | False] (default: False)
            True if head is to be counted in calculating the length.

          *head_width*: float or None (default: 3*width)
            total width of the full arrow head

          *head_length*: float or None (default: 1.5 * head_width)
            length of arrow head

          *shape*: ['full', 'left', 'right'] (default: 'full')
            draw the left-half, right-half, or full arrow

          *overhang*: float (default: 0)
            fraction that the arrow is swept back (0 overhang means
            triangular shape). Can be negative or greater than one.

          *head_starts_at_zero*: [True | False] (default: False)
            if True, the head starts being drawn at coordinate 0
            instead of ending at coordinate 0.

        Other valid kwargs (inherited from :class:`Patch`) are:
        %(Patch)s

        """
        if head_width is None:
            head_width = 20 * width
        if head_length is None:
            head_length = 1.5 * head_width

        distance = np.sqrt(dx ** 2 + dy ** 2)
        if length_includes_head:
            length = distance
        else:
            length = distance + head_length
        if not length:
            verts = []  # display nothing if empty
        else:
            # start by drawing horizontal arrow, point at (0,0)
            hw, hl, hs, lw = head_width, head_length, overhang, width
            left_half_arrow = np.array([
               [0.0, 0.0],                  # tip
                [-hl, -hw / 2.0],             # leftmost
                [-hl * (1 - hs), -lw / 2.0],  # meets stem
                [-length, -lw / 2.0],          # bottom left
                [-length, 0],
            ])
            #if we're not including the head, shift up by head length
            if not length_includes_head:
                left_half_arrow += [head_length, 0]
            #if the head starts at 0, shift up by another head length
            if head_starts_at_zero:
                left_half_arrow += [head_length / 2.0, 0]
            #figure out the shape, and complete accordingly
            if shape == 'left':
                coords = left_half_arrow
            else:
                right_half_arrow = left_half_arrow * [1, -1]
                if shape == 'right':
                    coords = right_half_arrow
                elif shape == 'full':
                    # The half-arrows contain the midpoint of the stem,
                    # which we can omit from the full arrow. Including it
                    # twice caused a problem with xpdf.
                    coords = np.concatenate([left_half_arrow[:-1],
                                             right_half_arrow[-2::-1]])
                else:
                    raise ValueError("Got unknown shape: %s" % shape)
            cx = float(dx) / distance
            sx = float(dy) / distance
            M = np.array([[cx, sx], [-sx, cx]])
            verts = np.dot(coords, M) + (x + dx, y + dy)

        Polygon.__init__(self, list(map(tuple, verts)), closed=True, **kwargs)

class MyArrow_Improvement(Polygon):
    """
    Like Arrow, but lets you set head width and head height independently.
    """
    def __str__(self):
        return "MyArrow()"

    def __init__(self, x, y, dx, dy, width=0.001, length_includes_head=False,
        head_width=None, head_length=None, shape='full', overhang=0,
        head_starts_at_zero=False, **kwargs):
        """
        Constructor arguments
          *width*: float (default: 0.001)
            width of full arrow tail

          *length_includes_head*: [True | False] (default: False)
            True if head is to be counted in calculating the length.

          *head_width*: float or None (default: 3*width)
            total width of the full arrow head

          *head_length*: float or None (default: 1.5 * head_width)
            length of arrow head

          *shape*: ['full', 'left', 'right'] (default: 'full')
            draw the left-half, right-half, or full arrow

          *overhang*: float (default: 0)
            fraction that the arrow is swept back (0 overhang means
            triangular shape). Can be negative or greater than one.

          *head_starts_at_zero*: [True | False] (default: False)
            if True, the head starts being drawn at coordinate 0
            instead of ending at coordinate 0.

        Other valid kwargs (inherited from :class:`Patch`) are:
        %(Patch)s

        """
        if head_width is None:
            head_width = 20 * width
        if head_length is None:
            head_length = 1.5 * head_width

        distance = n.sqrt(dx ** 2 + dy ** 2)
        if length_includes_head:
            length = distance
        else:
            length = distance + head_length
        if not length:
            verts = []  # display nothing if empty
        else:
            # start by drawing horizontal arrow, point at (0,0)
            hw, hl, hs, lw = head_width, head_length, overhang, width
            left_half_arrow = n.array([
               [0.0, 0.0],                  # cneter point
                [-hl, -hw / 2.0],             # the "feather", i.e. swept back point
                [-hl * (1 - hs), -lw / 2.0],  # the tip
                [-length, -lw / 2.0],          # outside corner of the tail
                [-length, 0],                   # symmetry corner of the tail (on y=0)
                ])
            # Now correct 
            #tip slope
            m = (left_half_arrow[2,1]-left_half_arrow[1,1])/(left_half_arrow[2,0]-left_half_arrow[1,0])
            # Seeking the tip to be on y = 0....
            x1,y1 = left_half_arrow[1]
            y2 = 0
            x2 = x1 + (y2-y1)/m
            left_half_arrow[2] = [x2,y2]
            #Move "origin" center point to intersect y=0
            #slope
            m = (left_half_arrow[0,1]-left_half_arrow[1,1])/(left_half_arrow[0,0]-left_half_arrow[1,0])
            # Seeking the tip to be on y = 0....
            x1,y1 = left_half_arrow[1]
            y2 = left_half_arrow[3,1]
            x2 = x1 + (y2-y1)/m
            left_half_arrow[0] = [x2,y2]
            #Now reorder
            left_half_arrow = left_half_arrow[[4,3,0,1,2]]
            
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
                    coords = n.vstack( (coords,coords[0]) )
                else:
                    raise ValueError("Got unknown shape: %s" % shape)
            cx = float(dx) / distance
            sx = float(dy) / distance
            M = n.array([[cx, sx], [-sx, cx]])
            verts = n.dot(coords, M) + (x + dx, y + dy)
            #codes = 
        Polygon.__init__(self, list(map(tuple, verts)), closed=True, **kwargs)

def myax_obsolete(fig,OH=2.9,HL=.012,HW=.018,TW=.003,AL=.14,PLW=0.0,rast=False):
    '''
    Converts a figure with single y-axis scale to K style.
    Requires figure hande or p.gcf()
    Other args are fancy arrow props.  
    OH(=2.9) is overhang
    HL(=.012) is head length.
    HW(=.018) is head width.
    TW(=.00p3) is tail width.
    AL(=.14) is arrow length.
    PLW(=0.0) is the patch kwarg linewidth
    rast(=False) whether to rasterise the patch.
     for arrow patch pertain to the y-axis arrow, and are scaled acoordingly for the x-axis.
     '''
    #Handle for figure axes
    ax = fig.gca()
    # Figure dimensions and axes bounds points in figure coordinates
    figwd,fight = fig.get_size_inches()
    bboxAX = fig.transFigure.inverted().transform(ax.bbox.get_points())
    x1,y1,x2,y2 = bboxAX.flatten()
    h2w = (y2-y1)*fight/(x2-x1)/figwd   #Ht to wt. conversion
    for i in ax.get_children():
        if type(i) is MyArrow:
            i.remove()
    #An inverse transform that will take me from display to data coordinates
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
            lab.set_size(30)

    # Now begin the arrow process
    # Get the bounding box for the x tick labels and convert it to axes coords, and set position accordingly
    bboxX = inv.transform( axX.get_ticklabel_extents(R)[0].get_points() )
    xlab.set_va('top')
    xlab.set_horizontalalignment('left')
    axX.set_label_coords( (3/4), n.min(bboxX[:,1]), transform = ax.transAxes)
    bboxLx = inv.transform(xlab.get_window_extent(R).get_points())
    left = n.min(bboxLx[:,0]) #Could cause problems! 
    OH1=OH*h2w
    HL1=HL*h2w
    HW1=HW*(1/h2w)
    TW1=TW*(1/h2w)
    AL1 = AL #(no *h2w because we want arrow length to scale with axis size)
    ar=MyArrow(left-4*HL1-AL1,n.mean(bboxLx[:,1]),AL1,0,overhang=OH,transform=ax.transAxes,
                    color='k',head_length=HL1,length_includes_head=True,
                    head_width=HW1,width=TW1,lw=PLW,rasterized=rast)
    ax.add_patch(ar)
    ar.set_clip_on(False)
    
    ########## Y label #############
    # Get bbox for ytick labels
    bboxY = inv.transform( axY.get_ticklabel_extents(R)[0].get_points() )
    # Rotate and postion ylabel, then translate it based on the bbox overlap
    ylab.set_verticalalignment('center')
    ylab.set_horizontalalignment('center')
    ylab.set_rotation(0)
    #xpos, ypos = n.min(bboxY[:,0]), (3/4) * n.abs(n.diff(bboxY[:,1])[0])
    xpos, ypos = n.min(bboxY[:,0]), (3/4) 
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes )
    bboxLy = inv.transform(ylab.get_window_extent(R).get_points())
    xpos = xpos - ( n.max(bboxLy[:,0]) - xpos )
    axY.set_label_coords( xpos , ypos , transform = ax.transAxes)
    bot = n.min(bboxLy[:,1])
    # Add arrow
    OH1=OH
    HL1=HL
    HW1=HW
    TW1=TW
    AL1=AL
    ar=MyArrow(xpos,bot-AL1-4*HL1,0,AL1,overhang=OH1,transform=ax.transAxes,
                    color='k',head_length=HL1,length_includes_head=True,
                    head_width=HW1,width=TW1,lw=PLW,rasterized=rast)
    ax.add_patch(ar)
    ar.set_clip_on(False)
    ar.set_clip_on(False)
        
    p.draw()
    None

def rightax_obsolete(fig,conversion,ylabel):
    """
    Adds a 2nd axis to right hand side.  Plots the label.  Returns handle for axes.
    """
    axL = fig.gca()
    axR = axL.twinx()
    y1, y2 = axL.get_ylim()
    axR.set_ylim(conversion(y1), conversion(y2))
    axR.set_ylabel(ylabel)
    axR.yaxis.set_label_coords(1.14,.7)
    #axR.yaxis.get_label().get_position()
    axR.yaxis.get_label().set_rotation(0)
    axR.yaxis.get_label().set_verticalalignment('bottom')
    axR.figure.canvas.draw()
    return axL, axR
