"""
These are convenience functions used by the particle tracking dynamics
code "Ldyn_".
"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.path as mpltPath
import pickle
import numpy as np

# utility function for making segment lists
def get_seg_list(X,N):
    sl = [X+str(n) for n in range(1,N+1)]
    return sl
    
# make a dict of handy segment lists
seg_dict = {'HoodCanalInner': ['H'+str(n)+'_s' for n in range(3,9)],
            'HoodCanal': get_seg_list('H',8),
            'SouthSound': get_seg_list('S',4),
            'Whidbey': get_seg_list('W',4),
            'PS': get_seg_list('A',3) + get_seg_list('M',6)
                    + get_seg_list('T',2) + get_seg_list('S',4)
                    + get_seg_list('H',8) + get_seg_list('W',4),
            'PSTrim': ['A2', 'A3'] + get_seg_list('M',6)
                    + get_seg_list('T',2) + get_seg_list('S',4)
                    + get_seg_list('H',8) + ['W1','W2'],
            'AISouth': ['A2', 'A3'],
            'PS_no_AI': get_seg_list('M',6)
                    + get_seg_list('T',2) + get_seg_list('S',4)
                    + get_seg_list('H',8) + get_seg_list('W',4),
            'SoG': get_seg_list('G',6),
            'Salish': get_seg_list('A',3) + get_seg_list('M',6)
                    + get_seg_list('T',2) + get_seg_list('S',4)
                    + get_seg_list('H',8) + get_seg_list('W',4)
                    + get_seg_list('G',6),
            }
            
def packer(x,y):
    # creates 2D array out of vectors, with columns x,y
    N = len(x)
    xy = np.concatenate((x.reshape(N,1), y.reshape(N,1)), axis=1)
    return xy
    
def get_imask(Ldir, seg_list, plon, plat, glon, glat):
    
    # load the TEF segment ji_dicts and volumes
    voldir = Ldir['LOo'] + 'tef2/volumes_cas6/'
    ji_dict = pickle.load(open(voldir + 'ji_dict.p', 'rb'))
    
    j_dict = {}; i_dict = {}
    for seg_name in seg_list:
        jj = []; ii = []
        ji_list_full = ji_dict[seg_name]
        for ji in ji_list_full:
            jj.append(ji[0])
            ii.append(ji[1])
        jjj = np.array(jj)
        iii = np.array(ii)
        j_dict[seg_name] = jjj
        i_dict[seg_name] = iii
    
    # now find the convex hull around the segments, and the points in that polygon
    path_dict = {}
    for sn in seg_list:
        xx = glon[j_dict[sn], i_dict[sn]]
        yy = glat[j_dict[sn], i_dict[sn]]
        nseg = len(xx)
        points = packer(xx,yy)
        hull = ConvexHull(points)
        # these xy, yh are the points that define the polygon of the hull
        xh = points[hull.vertices,0]
        yh = points[hull.vertices,1]
        xyh = packer(xh, yh)
        path_dict[sn] = mpltPath.Path(xyh)
    
    # and now find indices into plon and plat of winners
    xyp = packer(plon, plat)
    mask = np.ones_like(plon)==0 # all false to start
    for sn in seg_list:
        path = path_dict[sn]
        is_inside = path.contains_points(xyp)
        mask = mask | is_inside
    # and make a list of the winning particle numbers
    imask = [i for i, val in enumerate(mask) if val]
    # and imask will be in sorted order already
    
    return mask, imask

