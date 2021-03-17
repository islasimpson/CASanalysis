import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plotbox(fig, x1,x2,y1,y2, color='red'):
    """plot a colored box"""
    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    box = Rectangle((0.,0.),1,1,facecolor=color)
    ax.axis('off')
    return ax
