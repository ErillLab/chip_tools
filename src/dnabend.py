from dnacurve import Model
from utils import matrix_mult, matrix_add,transpose,pairs,wc
from math import sin,cos,pi

for model_name in "STRAIGHT AAWEDGE CALLADINE TRIFONOV DESANTIS REVERSED".split():
    exec("%s = Model.%s" %(model_name.lower(),model_name))

def functionify_model(model):
    """Accept a model defined in a dictionary and return a function
    from dinucleotides to properties"""
    def model_f(dinuc):
        print dinuc
        oligos = model["oligo"].split()
        if dinuc in oligos:
            i = oligos.index(dinuc)
            d = {prop:model[prop][i]
                 for prop in "twist tilt roll".split()}
        else:
            print "elsing"
            cunid = wc(dinuc)
            print cunid
            i = oligos.index(cunid)
            d = {prop:model[prop][i]
                 for prop in "twist roll".split()}
            print "d"
            d["tilt"] = -model["tilt"][i] # flip sign of tilt if reverse complementing
        d["rise"] = model["rise"]
        print "returning d"
        return d
    return model_f

def roll_matrix(roll):
    """Return a rotation matrix for roll- rotation in XZ plane"""
    phi = pi/180.0 * roll # convert degrees to rads
    return [[cos(phi),  0, sin(phi)],
            [0,         1,        0],
            [-sin(phi), 0, cos(phi)]]

def tilt_matrix(tilt):
    """Return a rotation matrix for tilt- rotation in YZ plane"""
    phi = pi/180.0 * tilt # convert degrees to rads
    return [[1,         0,        0],
            [0,  cos(phi), sin(phi)],
            [0, -sin(phi), cos(phi)]]

def twist_matrix(twist):
    """Return a rotation matrix for twist- rotation in XY plane"""
    phi = pi/180.0 * twist # convert degrees to rads
    return [[cos(phi),  sin(phi), 0],
            [-sin(phi), cos(phi), 0],
            [0,         0,        1]]

def main_example():
    sequence="""CGAAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAACGCG
AAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAACGCGAAAAAA
CGCGAAAAAACGCGAAAAAACG""".replace("\n","")
    lookup = functionify_model(aawedge)
    rise = 3.38 #Angstroms, from Gohlke
    b0 = transpose([[0,0,0]])
    v0 = transpose([[0,0,rise]])
    bs = [b0]
    vs = [v0]
    for pair in pairs(sequence):
        dinuc = "".join(pair)
        params = lookup(dinuc)
        ro, ti, tw = params["roll"],params["tilt"],params["twist"]
        last_b = bs[-1]
        last_v = vs[-1]
        new_v = reduce(matrix_mult,[roll_matrix(ro),
                                    tilt_matrix(ti),
                                    twist_matrix(tw),
                                    last_v])
        new_b = matrix_add(last_b,new_v)
        bs.append(new_b)
        vs.append(new_v)
    return bs

def distance(u,v):
    [[x1],[y1],[z1]],[[x2],[y2],[z2]] = u,v
    return sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def distances(bs):
    return map(uncurry(distance),pairs(bs))

def test_plot():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    xs,ys,zs = map(concat,transpose(main_example()))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs)
    plt.show()
    


                       
print "loaded"
