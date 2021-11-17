import numpy as np
from pymatgen.core import Structure

class PhonopySumoProjections:
    ''' grab the phonons using the sumo interface '''
    def __init__(self,
                 FORCE_SETS='FORCE_SETS',
                 unitcell = 'POSCAR',
                 dim = None):
        
        self.forces = FORCE_SETS
        self.unitcell = unitcell
        self.structure = Structure.from_file(unitcell)
        self.dim = [[dim[0],0,0],
                    [0,dim[1],0],
                    [0,0,dim[2]]]
        
    def _get_sumo_phonopy_object(self,**kwargs): # sumo takes the edge out of getting the data together
        from sumo.cli.phonon_bandplot import _bs_from_filename as sumo_bs_from_filename
        bs, phon = sumo_bs_from_filename(filename=self.forces, 
                                         poscar = self.unitcell, 
                                         dim = self.dim,
                                         symprec = 0.01,
                                         spg = None, #self.structure.get_space_group_info()[1]
                                         kpt_list = None,
                                         labels = None,
                                         primitive_axis = None,
                                         units = 'THz',
                                         born = False,
                                         mode = 'bradcrack',
                                         eigenvectors = True,
                                         line_density = 100)
        
        return(bs)
    
    def _get_elemental_phonon_weights(self,element, bs,**kwargs):
        ''' takes one element at a time '''
        from pymatgen.io.phonopy import eigvec_to_eigdispl
        from pymatgen.core.periodic_table import Element
        
        nq = bs.nb_qpoints #number of qpoints
        nb = bs.nb_bands #number of bands
        na = bs.structure.num_sites # number of atoms/sites
        vs = [i for i,x in enumerate(bs.structure.species) if x == Element(element)] # valid sites 
        
        weights = []
        for b in range(nb):
            bnds = []
            for q in range(nq):
                qps = []
                for a in range(na):
                    if a in vs:
                        #print(q,b,a)
                        ats = np.linalg.norm(eigvec_to_eigdispl(v=bs.eigendisplacements[b][q][a],
                                                                q=bs.qpoints[q].frac_coords,
                                                                frac_coords=bs.structure[a].frac_coords,
                                                                mass=Element(element).atomic_mass))
                        
                        qps.append(ats)
                bnds.append(qps)
            weights.append(bnds)

        return(np.asarray(weights))
    
    def create_plot_scatter(self,element,ax,cmap='Blues',alpha=0.3,size=10,**kwargs):
        from sumo.plotting import sumo_base_style
        import matplotlib.pyplot as plt
        plt.style.use(sumo_base_style) # this can probably be done in a better way 
        from matplotlib.colors import Normalize
        
        bs = self._get_sumo_phonopy_object()
        weights = self._get_elemental_phonon_weights(element=element,bs=bs)
        
        ax = ax 
        
        nq = bs.nb_qpoints
        nb = bs.nb_bands
        
        ax.set_xlim(0,nq)
        norm = Normalize(vmin=np.min(weights),vmax=np.max(weights))
        q = np.arange(bs.nb_qpoints)
        for i in np.arange(nb):
            y = [y[i] for y in bs.bands.T]
            w = weights[i]
            ax.scatter(q,y,c=w,s=size,cmap=cmap,alpha=alpha,edgecolors='None',norm=norm,**kwargs)
            
        ax.set_ylabel('Frequency (THz)')

        tick_points = []
        tick_labels = []
        for i,k in enumerate(bs.qpoints):
            if not k.label == None:
                if not i == 0:
                    if bs.qpoints[i-1].label == None:
                        tick_points.append(i)
                        if k.label == '\\Gamma':
                            tick_labels.append('$\\Gamma$')
                        else:
                            tick_labels.append(k.label)
                else:
                    if k.label == '\\Gamma':
                        tick_labels.append('$\\Gamma$')
                        tick_points.append(i)
                    else:
                        tick_labels.append(k.label)

        ax.set_xticks(tick_points)
        ax.set_xticklabels(tick_labels)
        [ax.axvline(i,color='k') for i in tick_points]
        return(ax)
