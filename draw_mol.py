#!/usr/bin/env python
import rdkit.Chem.Draw as Draw
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem as Chem
import matplotlib.cm as cm

if __name__ == "__main__":
    m = Chem.MolFromSmiles("CCC")
    fig=Draw.MolToMPL(m)
    x,y,z=Draw.calcAtomGaussians(m,0.03,step=0.01,weights=(40,1,3))
    fig.axes[0].imshow(z,cmap=cm.Oranges,interpolation='bilinear',origin='lower',extent=(0,1,0,1))
    fig.axes[0].contour(x,y,z,20,colors='k',alpha=0.5)
    #fig.show()
    fig.savefig('coumlogps.colored.png',bbox_inches='tight')
