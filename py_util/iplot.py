'''
ipython/pylab helpers for plotting

'''
import pylab,numpy,matplotlib,os,sys,re
from copy import copy

histcolors = {'b':'blue',
              'g':'green',
              'r':'red',
              'c':'cyan',
              'm':'magenta',
              'y':'yellow',
              'k':'black',
              'w':'white',}

dig = range(9)+['A','B','C','D','E','F']

spectrum = ['#FF%s%s00' % (i,j) for i in dig for j in dig] + ['#%s%sFF00' % (i,j) for i in dig[::-1] for j in dig[::-1]] + ['#00FF%s%s' % (i,j) for i in dig for j in dig] +['#00%s%sFF' % (i,j) for i in dig[::-1] for j in dig[::-1]] + ['#%s%s00FF' % (i,j) for i in dig for j in dig]

spectrum_noUV = ['#FF%s%s00' % (i,j) for i in dig for j in dig] + ['#%s%sFF00' % (i,j) for i in dig[::-1] for j in dig[::-1]] + ['#00FF%s%s' % (i,j) for i in dig for j in dig] +['#00%s%sFF' % (i,j) for i in dig[::-1] for j in dig[::-1]]

def to_greyscale(hex_str):
    '''takes a hex color string ("#FF0000") and converts to grayscale'''
    col = int(hex_str.replace('#','0x'),16)
    return hex((((((((col >> 16) & 0xff)*76) + (((col >> 8) & 0xff)*150) +((col & 0xff)*29)) >> 8)) << 16) |(((((((col >> 16) & 0xff)*76) + (((col >> 8) & 0xff)*150) + ((col & 0xff)*29)) >> 8)) << 8) | ((((((col >> 16) & 0xff)*76) + (((col >> 8) & 0xff)*150) + ((col & 0xff)*29)) >> 8))).replace('0x','#')

def subspectrum(numcol,reversed=False,noUV=False):
    '''given a number of colors desired (<= len spectrum)
    returns a list of hex codes for colors'''
    #pad spectrum if undersized
    if noUV:
        spec = copy(spectrum_noUV)
    else:
        spec = copy(spectrum)
    if reversed:
        spec.reverse()
    while len(spec) < numcol:
        spec = list(reduce(lambda x,y:x+y,zip(spec,spec)))
    try:
        #step = len(spec)/numcol
        #return spec[::step][:numcol]
        return [spec[i] for i in map(int,numpy.linspace(0,len(spec)-1,numcol))]
    except:
        return []

def subspec_enum(iterable):
    return zip(subspectrum(len(iterable)),iterable)


def violin_plot(ax,data,pos, bp=False):
    '''
    create violin plots on an axis
    '''
    from scipy.stats import gaussian_kde
    from numpy import arange
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        #k.set_bandwidth(bw_method=k.factor / 2.)
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=pos,vert=1)


def draw_line_hist(data,bins=20,linewidth=3,color=None,normed=False,fig=1):
    '''draws an outline histogram, returns values,bins
    '''

    pylab.figure(fig)
    values,bins = pylab.histogram(data,bins=bins,new=True,normed=normed)
    if color:
        pylab.plot(bins[:-1],values,color,drawstyle='steps-post',linewidth=linewidth)
    else:
        pylab.plot(bins[:-1],values,drawstyle='steps-post',linewidth=linewidth)        
              
    return values,bins
    

def draw_hist_cdf(data,fig=None,nbins=None,subpl=None,nohist=False,**figargs):
    '''input data is a list of dicts with keys: data,histcolor,plotline
    '''

    bins = None

    if fig: pylab.figure(fig,**figargs)
    if subpl: pylab.subplot(subpl)
    
    results = []
    for d in data:
        if bins is not None:
            n,bins,patches=pylab.hist(d['data'],bins=bins,normed=1,fc=histcolors[d['plotline'][0]])
        elif nbins is not None:
            n,bins,patches=pylab.hist(d['data'],bins=nbins,normed=1,fc=histcolors[d['plotline'][0]])
        else:
            n,bins,patches=pylab.hist(d['data'],normed=1,fc=histcolors[d['plotline'][0]])
        results.append(n)

    if nohist:
        pylab.cla()
    else:
        pylab.twinx()

    for i,d in enumerate(data):
        y = pylab.plot(bins,numpy.cumsum(results[i])/sum(results[i]), d['plotline'], linewidth=2)

def windowed_mapping_from_gff(gffs,seqlens=None,winsize=500,winstep=250,scorefunct=lambda li: reduce(lambda x,y: x+y,li)):
    '''generates a windowed paired mapping gff from paired mapping gff

    if seqlens is not given will map until the last end in each
    score for each window is scorefunct([all_scores_in_window])    
    returns mapping gff
    '''
    
    if seqlens is None:
        seqlens = (max([r['end'] for r in gffs[0]]),
                   max([r['end'] for r in gffs[1]]),
                   )

    new_gffs = (GFF.File(),GFF.File())
    n = 1
    for q in xrange(1,seqlens[0],winstep):
        for s in xrange(1,seqlens[1],winstep):

            Qstart,Qend = q, q+winsize
            Sstart,Send = s, s+winsize
            print >>sys.stderr,Qstart,Qend,'->',Sstart,Send
            scores = [r['score'] for i,r in enumerate(gffs[0]) if (Qstart <= r['start'] 
                                                                   and Qend >= r['end']
                                                                   and Sstart <= gffs[1][i]['start']
                                                                   and Send >=gffs[1][i]['end']
                                                                   )]
            if scores:
                new_score = scorefunct(scores)
                Qreg = GFF.Region(fields={'seqid':gffs[0][0]['seqid'],'start':Qstart,'end':Qend,'score':new_score,'attribute_ID':str(n)})
                Sreg = GFF.Region(fields={'seqid':gffs[1][0]['seqid'],'start':Sstart,'end':Send,'score':new_score,'attribute_ID':str(n)})
                print >>sys.stderr,'adding',Qreg
                print >>sys.stderr,'adding',Sreg
                new_gffs[0].append(Qreg)
                new_gffs[1].append(Sreg)
                n += 1

    return new_gffs

def gff_to_mapping_data(gffs,ranks=(2,1),min_opacity=0.2,min_topscore=None):
    '''converts paired GFF lists into maps suitable for plotting

    ranks are currently in descending order
    '''
    scores = numpy.array([r['score'] for r in gffs[0]])

    newmax = max(scores.max(),min_topscore)
    shiftdown = scores.min() - (min_opacity*newmax)
    scores = (scores - shiftdown) / newmax

    data = []
    for i,reg in enumerate(gffs[1]):
        if reg['strand'] != gffs[0][i]['strand']:
            data.append((
                (gffs[0][i]['start'],gffs[0][i]['end'],reg['start'],reg['end'],ranks[0],ranks[1]),
                'red',
                scores[i]))
        else:
            data.append((
                (gffs[0][i]['start'],gffs[0][i]['end'],reg['start'],reg['end'],ranks[0],ranks[1]),
                'blue',
                scores[i]))
    return data

#from Aln import dotplot
#dotplot moved here, but better to use lastz -self, visualize w/ R

def dotplot(seq1=None,seq2=None,matrices=None,
            filebase=None,filetype='.png',
            window=20, cut=0.6, rev=False,normed=True,pylab_colors='gray',extend=False):

    colorset = getattr(pylab,pylab_colors)
    colorset()

    if filebase:
        outdir = os.path.dirname(filebase)
        if not os.path.exists(outdir): os.makedirs(outdir)

        if seq1 and seq2: filebase += '_%s-%0.2f' % (window,cut)
        if rev:
            filebase += '_rev'
        if normed:
            filebase += '_normed'
        if extend:
            filebase += '_extend'


    if seq1 and seq2:
        if os.path.exists(filebase+'-f.pkl') and os.path.exists(filebase+'-r.pkl'):
            f = cPickle.load(open(filebase+'-f.pkl'))
            r = cPickle.load(open(filebase+'-r.pkl'))
        else:
            (f,r) = bothstrands_id_matrix(seq1,seq2,rev=rev,cut=cut,window=window,normed=normed)
            cPickle.dump(f,open(filebase+'-f.pkl','wb'))
            cPickle.dump(r,open(filebase+'-r.pkl','wb'))

    elif matrices:
        if isinstance(matrices[0],(str,unicode)):
            f = cPickle.load(open(matrices[0]))
        else:
            f = matrices[0]
        if len(matrices) > 1:
            if isinstance(matrices[1],(str,unicode)):
                r = cPickle.load(open(matrices[1]))
            else:
                r = matrices[1]
        else:
            if rev:
                r = numpy.ones([len(f),len(f[0])],float)
            else:
                r = numpy.zeros([len(f),len(f[0])],float)

    else:
        raise ValueError, 'either seq1 AND seq2 must be set, or at least one matrix or matrix file must be supplied'


    if extend:
        #extend fwd hits:
        if rev:
            f_mask = numpy.ones([len(f),len(f[0])],float)
        else:
            f_mask = numpy.zeros([len(f),len(f[0])],float)
        for i in range(len(f)-window*2):
            for j in range(len(f[i])-window*2):
                for cnt in range(window):
                    if rev:
                        if f_mask[i+cnt][j+cnt] > f[i][j]: f_mask[i+cnt][j+cnt] = f[i][j]
                    else:
                        if f_mask[i+cnt][j+cnt] < f[i][j]: f_mask[i+cnt][j+cnt] = f[i][j]
        f = f+f_mask
        #extend rev hits
        if rev:
            r_mask = numpy.ones([len(r),len(r[0])],float)
        else:
            r_mask = numpy.zeros([len(r),len(r[0])],float)
        for i in range(len(r)-window):
            for j in range(len(r[i])-1,window,-1):
                for cnt in range(window):
                    if rev:
                        if r_mask[i+cnt][j-cnt] > r[i][j]: r_mask[i+cnt][j-cnt] = r[i][j]
                    else:
                        if r_mask[i+cnt][j-cnt] < r[i][j]: r_mask[i+cnt][j-cnt] = r[i][j]
        r = r+r_mask
        

    fig = pylab.matshow(f+r)
    if filebase:
        pylab.title(os.path.basename(filebase))
        pylab.savefig(filebase+'_'+pylab_colors+filetype)
        
    return (f,r,fig)


import GFF
from copy import deepcopy
def merge_mapping_overlaps(data):
    '''given a mapping list of tuples data, should merge any two mappings on the same strand, that overlap in both species

    '''

    newdata = deepcopy(data)
    for i,(this_coords,this_color,this_score) in enumerate(newdata):
        for j,(other_coords,other_color,other_score) in enumerate(newdata):
            if ((this_color == other_color) and 
                ((other_coords[0] < this_coords[1] < other_coords[1]) or 
                 (this_coords[0] < other_coords[1] < this_coords[1])) and 
                ((other_coords[2] < this_coords[3] < other_coords[3]) or 
                 (this_coords[2] < other_coords[3] < this_coords[3]))):
                new_score = max(this_score,other_score)
                new_coords = (min(this_coords[0],other_coords[0]),
                              max(this_coords[1],other_coords[1]),
                              min(this_coords[2],other_coords[2]),
                              max(this_coords[3],other_coords[3]),) + this_coords[4:]
                newdata.append((new_coords,this_color,new_score))
                del newdata[j]
                break
        del newdata[i]
    return newdata
            

def draw_mapping_plot(data,seqbounds,fig=1,subpl=111,opacity_cutoff=0,renorm=False,bar_thickness=0.05,filename=None,ylabel=None,clear_plot=False,merge_overlaps=False,**figargs):
    '''input: 
    data is a list or tuples -- tuple(s1,e1,s2,e2),color,alpha: 
    [((start1,end1,start2,end2,rank1,rank2),color,opacity,), ... ]

    seqbounds is a 2-list of 2-tuples 
    [(seq1 start, seq1 stop),(seq2 start, seq2 stop)]
    OR a tuple of 2 numbers
    (seq1 stop,seq2 stop)
    
    returns figure object
    '''

    if isinstance(seqbounds,tuple):
        seqbounds = [(1,seqbounds[0]),(1,seqbounds[1])]

    print >> sys.stderr,'seqbounds',seqbounds 
    if merge_overlaps:
        print >> sys.stderr, 'merged',len(data),
        data = merge_mapping_overlaps(data)
        print >> sys.stderr, 'to',len(data)

    figo = pylab.figure(fig,**figargs)
    if clear_plot:
        figo.clf()
        figo = pylab.figure(fig,**figargs)

    ax = figo.add_subplot(subpl)
    ax.set_yticks([])
    if ylabel: ax.set_ylabel(ylabel)
    
    #patch this to use offsets passed in

    offsets = [0,0]
    [(s1_start,s1_end),(s2_start,s2_end)] = seqbounds
    s1_len = s1_end-s1_start
    s2_len = s2_end-s2_start
    if s1_len > s2_len:
        offsets[1] = s1_len/2 - s2_len/2
    else:
        offsets[0] = s2_len/2 - s1_len/2

    print >> sys.stderr,'offsets',offsets

    (rank1,rank2) = (None,)*2

    print >> sys.stderr,'data',len(data),'above cutoff',opacity_cutoff,
    data = [d for d in data if d[2] > opacity_cutoff]
    print >> sys.stderr,len(data)

    if renorm:
        scores = numpy.array([d[2] for d in data])
        max_val = scores.max()
        min_val = scores.min()

    for coords,color,opacity in data:
        (rank1,rank2) = coords[4:6]

        if renorm: 
            opacity = (opacity - min_val) / max_val

        top = rank1-bar_thickness
        bot = rank2+bar_thickness
        coords = (coords[0]+offsets[0],
                  coords[1]+offsets[0],
                  coords[3]+offsets[1],
                  coords[2]+offsets[1])
        points=zip(coords[:4],((top,) * 2 + (bot,) * 2))
        ax.add_patch(matplotlib.patches.Polygon(points,alpha=opacity,fc=color,ec=color))

    if bar_thickness:
        print >> sys.stderr, 'setting up bars:',s1_end,rank1,s2_end,rank2

        points = [(offsets[0],rank1),
                  (offsets[0]+s1_end,rank1),
                  (offsets[0]+s1_end,rank1-bar_thickness),
                  (offsets[0],rank1-bar_thickness)]
        ax.add_patch(matplotlib.patches.Polygon(points,fc='w',ec='k'))
        
        points = [(offsets[1],rank2+bar_thickness),
                  (offsets[1]+s2_end,rank2+bar_thickness),
                  (offsets[1]+s2_end,rank2),
                  (offsets[1],rank2)]
        ax.add_patch(matplotlib.patches.Polygon(points,fc='w',ec='k'))
    
    if filename:

        outdir = os.path.dirname(filename)
        if not os.path.exists(outdir): os.makedirs(outdir)
        pylab.title(os.path.basename(filename))

        pylab.plot()
        figo.savefig(filename)
    else:
        pylab.plot()
    return figo

def mapping_from_array(id_array,color,window,rev=False):
    d = []
    
    id_array = id_array - id_array.min()
    id_array = id_array / id_array.max()
    if rev:
        id_array = 1-id_array
        
    for s1, s2, score in ((i,j,id_array[i][j]) for i,vi in enumerate(id_array) for j,vj in enumerate(id_array[i])):
        if score:
            d.append(
                ((s1,s1+window,s2,s2+window,2,1),color,score) )
    return d

   
def mapping_from_id_arrays(id_arrays,window,fig=1,opacity_cutoff=0,renorm=False,bar_thickness=0.05,filename=None,clear_plot=True,rev=False,merge_overlaps=False,**figargs):


    #hack to clear figure
    figo = pylab.figure(fig,**figargs)
    if clear_plot:
        figo.clf()
        figo = pylab.figure(fig,**figargs)

    for i,id_array in enumerate(id_arrays):
        subpl = int(str(len(id_arrays))+'1'+str(i+1))

        data = []
        if isinstance(id_array,tuple):
            if len(id_array) > 2:
                ylabel = id_array[2]
            else:
                ylabel = None

            f_data = []
            r_data = []
            if id_array[0].max()-id_array[0].min(): f_data = mapping_from_array(id_array[0],'blue',rev=rev,window=window)
            if id_array[0].max()-id_array[1].min(): r_data = mapping_from_array(id_array[1],'red',rev=rev,window=window)
            data = r_data + f_data
            
            seqbounds = [(0,len(id_array[0])+window),(0,len(id_array[0][0])+window)]
        else:
            ylabel=None
            if id_array[0].max()-id_array[0].min(): data.extend(mapping_from_array(id_array,'black',rev=rev,window=window))
            seqbounds = [(0,len(id_array)+window),(0,len(id_array[0])+window)]
            

        fig_drawn = draw_mapping_plot(data,seqbounds,opacity_cutoff=opacity_cutoff,renorm=renorm,subpl=subpl,fig=fig,ylabel=ylabel,bar_thickness=bar_thickness,merge_overlaps=merge_overlaps)

    if filename:
        pylab.subplot(int(str(len(id_arrays))+'11'))
        pylab.title(os.path.basename(filename))
        outdir = os.path.dirname(filename)
        if not os.path.exists(outdir): os.makedirs(outdir)

        
        figo.savefig(filename)
    return fig_drawn
import GFF, Util
def draw_bs_plot(sites,sp_order,site_styles,seq_lens,offsets=None,maxheight=0.8,minheight=0.4,
                 fig=1,subpl=111,clear_plot=True,filename=None,**figargs):

    by_factor = dict(zip(set([r['source'] for r in sites]),[GFF.File() for i in set([r['source'] for r in sites])]))

    for r in sites:
        cut = site_styles[r['source']]['cut']
        if r['score'] < cut and r['seqid'] in sp_order:
            by_factor[r['source']].append(r)


    print by_factor

    for k,v in by_factor.items():
        normscores = Util.normalize([r['score'] for r in v],minheight,maxheight,to_abs=1)
        for i,vn in enumerate(normscores):
            by_factor[k][i]['score'] = vn

    sites_to_plot = []
    for f in by_factor.values():
        sites_to_plot.extend(f)
            
    figo = pylab.figure(fig,**figargs)
    if clear_plot:
        figo.clf()
        figo = pylab.figure(fig,**figargs)

    ax = figo.add_subplot(subpl)
    ax.set_yticks([])

    #calc offsets, draw lines
    if offsets is None:
        offsets = [None]*(len(sp_order)+1)
        midpt = max([v for k,v in seq_lens.items() if k in sp_order])/2
    for i,sp in enumerate(sp_order):
        rank = len(sp_order) - i
        if offsets[rank] is None:
            off = midpt - seq_lens[sp]/2
            offsets[rank] = off
            print off,rank,seq_lens[sp]+off,rank
        ax.text(5,rank,sp)
        ax.add_line(matplotlib.lines.Line2D((offsets[rank],seq_lens[sp]+offsets[rank]),(rank,rank),color='k',alpha=0.25,lw=5))
        

    for site in sites_to_plot:
        fc = site_styles[site['source']]['color']
        ec = fc
        rank = len(sp_order) - sp_order.index(site['seqid'])
        ax.add_patch(matplotlib.patches.Ellipse( (site['start']+offsets[rank],rank),
                                                   len(site),
                                                   site['score'],
                                                   fc=fc,ec=ec,alpha=site['score'] )
                     )
    
    if filename:
        ax.autoscale_view()
        figo.savefig(filename)
    else:
        pylab.plot()

def draw_gff(gff,offset=0,rank=1,height=0.8,alpha_range=None,arrowhead=None,color='black',fig=1,subpl=111,clear_plot=False,**figargs):
    '''currently deprecated--use add_gff_patches until this can be re-written to do so!
    '''
    figo = pylab.figure(fig,**figargs)
    if clear_plot:
        figo.clf()
        figo = pylab.figure(fig,**figargs)
    
    ax = figo.add_subplot(subpl)
    ax.set_yticks([])
    
    y = rank - (height / 2)

    if alpha_range:
        normscores = Util.normalize([r['score'] for r in gff],alpha_range[0],alpha_range[1],to_abs=1)
        for i,r in enumerate(gff): #add strands, make funct
            if arrowhead:
                points = ( (r['start'],rank+(height/2)),
                           ((len(r)*(1-arrowhead))+r['start'],rank+(height/2)),
                           (r['end'],rank),
                           ((len(r)*(1-arrowhead))+r['start'],rank-(height/2)),
                           (r['start'],rank-(height/2))
                           )
            else:
                points = ( (r['start'],rank+(height/2)),
                           (r['end'],rank+(height/2)),
                           (r['end'],rank-(height/2)),
                           (r['start'],rank-(height/2))
                           )
            ax.add_patch(
                matplotlib.patches.Polygon( points, fc=color,ec=color,alpha=normscores[i] ))
    else:
        for r in gff:
            if arrowhead:
                points = ( (r['start'],rank+(height/2)),
                           ((len(r)*(1-arrowhead))+r['start'],rank+(height/2)),
                           (r['end'],rank),
                           ((len(r)*(1-arrowhead))+r['start'],rank-(height/2)),
                           (r['start'],rank-(height/2))
                           )
            else:
                points = ( (r['start'],rank+(height/2)),
                           (r['end'],rank+(height/2)),
                           (r['end'],rank-(height/2)),
                           (r['start'],rank-(height/2))
                           )
            print points
            ax.add_patch(
                matplotlib.patches.Polygon( points, fc=color,ec=color ))
            
    
def add_filled_plot(ax_obj,yvals,xvals=None,color='black',offset=0,rank=0,baseline=None,maxpoints=32000):
    '''adds a polygon of a filled plot (originally intended for phastCons scores and the like)
    to an axes instance (e.g. fig.axes[0])

    rank is the y value of 'baseline' on the plot (not the midpoint, as is the behavior of other poly generators)
    '''
    if len(yvals) <= maxpoints:
        yvals = [y+rank for y in yvals]
        
        if baseline is None:
            baseline = min(yvals)

        if xvals is None:
            xvals = range(offset,len(yvals)+offset)
            
        points = zip(xvals,yvals) + [(xvals[-1],baseline)] + [(offset,baseline)]
        poly = matplotlib.patches.Polygon(points,fc=color,ec=color)
        ax_obj.add_patch(poly)
        return [poly]
    else:
        split = len(yvals)/2
        
        if xvals:
            return add_filled_plot(ax_obj,yvals[:split+1],xvals[:split+1],color,offset,rank,baseline,maxpoints) + \
                add_filled_plot(ax_obj,yvals[split:],xvals[split:],color,offset=offset+split,rank=rank,baseline=baseline,maxpoints=maxpoints)
        else:
            return add_filled_plot(ax_obj,yvals[:split+1],xvals,color,offset,rank,baseline,maxpoints) + \
                add_filled_plot(ax_obj,yvals[split:],xvals,color,offset=offset+split,rank=rank,baseline=baseline,maxpoints=maxpoints)


def add_gff_patches(ax_obj,gff,offset=0,rank=1,height=0.8,alpha_range=None,arrowhead=None,gene_exon_types=None,color='black'):
    
    if gene_exon_types is not None and len(gene_exon_types) == 2:
        polys = []
        polys.extend(
            add_gff_patches(ax_obj,[r for r in gff if r['type'] == gene_exon_types[0]],offset,rank,height*0.1,alpha_range,arrowhead=arrowhead,color=color)
            )
        polys.extend(
            add_gff_patches(ax_obj,[r for r in gff if r['type'] == gene_exon_types[1]],offset,rank,height,alpha_range,arrowhead=arrowhead,color=color)
            )
        return polys
        
        

    def get_points(region,arrowhead,offset,rank):
        if arrowhead:
            if region['strand'] == '-':
                points = ( (r['end']+offset,rank+(height/2)),
                           ((len(r)*(arrowhead))+r['start']+offset,rank+(height/2)),
                           (r['start']+offset,rank),
                           ((len(r)*(arrowhead))+r['start']+offset,rank-(height/2)),
                           (r['end']+offset,rank-(height/2))
                           )

            else:
                points = ( (r['start']+offset,rank+(height/2)),
                           ((len(r)*(1-arrowhead))+r['start']+offset,rank+(height/2)),
                           (r['end']+offset,rank),
                           ((len(r)*(1-arrowhead))+r['start']+offset,rank-(height/2)),
                           (r['start']+offset,rank-(height/2))
                           )

        else:
            points = ( (r['start']+offset,rank+(height/2)),
                       (r['end']+offset,rank+(height/2)),
                       (r['end']+offset,rank-(height/2)),
                       (r['start']+offset,rank-(height/2))
                       )
        return points
        
    y = rank - (height / 2)
    polys = []
    if alpha_range:
        normscores = Util.normalize([r['score'] for r in gff],alpha_range[0],alpha_range[1],to_abs=1)
        for i,r in enumerate(gff): #add strands, make funct
            points = get_points(r,arrowhead,offset,rank)
            polys.append(matplotlib.patches.Polygon( points, fc=color,ec=color,alpha=normscores[i] ) )

    else:
        for r in gff:
            points = get_points(r,arrowhead,offset,rank)
            polys.append( matplotlib.patches.Polygon( points, fc=color,ec=color ) )

    for p in polys: ax_obj.add_patch(p)
    return polys

def add_splitbox(ax_obj,box_origin,ltcol,rtcol,slope=1,box_size=1):
	'''given an axes object, an origin (lower left corner of splitbox) two colors (left and right) and 
	'slope' which can be either -1 or 1, and will result in a splitline either:
	positive (bottom left to top right) or negative (top left to bottom right)
	
	'''
	lastlt = box_size*0.5 + slope*0.5
	lastrt = box_size*0.5 - slope*0.5
	ox,oy = box_origin
	#lt = [box_origin, (ox,oy+box_size), (ox+box_size,lastlt)]
	#rt = [(ox+box_size,oy+box_size), (ox+box_size,oy), (ox,lastrt)]
	lt = [box_origin, (ox,oy+box_size), (ox+box_size,oy+box_size)]
	rt = [box_origin, (ox+box_size,oy), (ox+box_size,oy+box_size)]
	
	ax_obj.add_patch(matplotlib.patches.Polygon(lt,fc=ltcol))
	ax_obj.add_patch(matplotlib.patches.Polygon(rt,fc=rtcol))
	
def load_cie_funcs(filename):
	cie = {}
	for l in open(filename):
		match = re.search(r'([\d\.]+),\s+(\-?[\d\.eE\+-]+),\s+(\-?[\d\.eE\+-]+),\s+(\-?[\d\.eE\+-]+)',l)
		if match:
			f = [float(i) for i in match.groups()]
			cie[int(f[0])] = f[1:]
	return cie
	
def calc_cie_XYZ(cie_funcs,intensities):
	fns = ['X','Y','Z']
	#cie_xyz = {}.fromkeys(fns,0)
	cie_xyz = [0.0,0.0,0.0]
	for l,i in intensities:
		for it,fn in enumerate(fns):
			cie_xyz[it] += cie_funcs.get(int(l),(0.0,0.0,0.0))[it]*i
	cie_xyz = [max(i,0) for i in cie_xyz]
	tot = sum(cie_xyz)
	for k,v in enumerate(cie_xyz):
		cie_xyz[k] = v/tot
	return cie_xyz

## {{{ Recipe 412982 (r1): Use PIL to make a "contact sheet" montage of images 
def make_contact_sheet(fnames,(ncols,nrows),(photow,photoh),
                       (marl,mart,marr,marb),
                       padding):
    """\
    Make a contact sheet from a group of filenames:

    fnames       A list of names of the image files
    
    ncols        Number of columns in the contact sheet
    nrows        Number of rows in the contact sheet
    photow       The width of the photo thumbs in pixels
    photoh       The height of the photo thumbs in pixels

    marl         The left margin in pixels
    mart         The top margin in pixels
    marr         The right margin in pixels
    marl         The left margin in pixels

    padding      The padding between images in pixels

    returns a PIL image object.
    """

    # Read in all images and resize appropriately
    imgs = [Image.open(fn).resize((photow,photoh)) for fn in fnames]

    # Calculate the size of the output image, based on the
    #  photo thumb sizes, margins, and padding
    marw = marl+marr
    marh = mart+ marb

    padw = (ncols-1)*padding
    padh = (nrows-1)*padding
    isize = (ncols*photow+marw+padw,nrows*photoh+marh+padh)

    # Create the new image. The background doesn't have to be white
    white = (255,255,255)
    inew = Image.new('RGB',isize,white)

    # Insert each thumb:
    for irow in range(nrows):
        for icol in range(ncols):
            left = marl + icol*(photow+padding)
            right = left + photow
            upper = mart + irow*(photoh+padding)
            lower = upper + photoh
            bbox = (left,upper,right,lower)
            try:
                img = imgs.pop(0)
            except:
                break
            inew.paste(img,bbox)
    return inew
## End of recipe 412982 }}}
