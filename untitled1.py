# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 14:35:21 2015

@author: jordan
"""
import matplotlib.pyplot as plt
import CreateCells as cc
import inputexplorere as ie

#def main():
ab = cc.Create_Cell('ABfile')

def model(*args):
    for i, x in enumerate(args):

        p = labels[i].partition('_')

        ab._optimize[p[0]][p[2]] = x
    ab.ga_change_gbar(ab._optimize)
    return ab.run(tstop=10000) 

sliders=[]
labels = []
for key, val in ab._optimize.items():
    for key2, val2 in val.items():
        sliders.append({'label':(key+'_'+key2),
        'valmin': val2[0],
        'valmax': val2[1],
        'valinit': ab._initials[key][key2]})
        labels.append(key+'_'+key2)

fig, ax = plt.subplots(1)
init = [x['valinit'] for x in [s for s in sliders]]
ax.plot(model(*init)['time'],model(*init)['vs'])
def dynamicplot(*args):
    ax.clear()
    ax.plot(model(*args)['time'],model(*args)['vs'])
    fig.canvas.draw()

ie.inputExplorer(dynamicplot,sliders)

#if __name__== '__main__':
#    main()