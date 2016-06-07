#!/usr/bin/python
import numpy as np
import math
import sys
import os
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath']  
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica', 'Arial', 'Verdana', 'sans-serif']})
###serif_fonts = [ 'New Century Schoolbook', 'Times New Roman', 'Times', 'Palatino', 'serif']

class ImplementationError(Exception):
  pass

class VaspFormatError(Exception):
  pass

class DOSPlotter:
  """
  The DOS Plotter class.
  """
  def __init__(self):
    self.spin = False
    self.bins = None
    self.tdict = None
    self.pdict = None
    self.natoms = None
    self._natoms = None
    self.elems = None
    self.projs = None

  def read_tdos(self, fdos='DOSCAR'):
    """
    Reads the total density of states (DOS) data from the specified DOSCAR file,
    and stores the array of energies and DOS in a dictionary tdict.
    """
    if not os.path.exists(fdos):
      print "The file you specified does not exist."
      sys.exit(0)

    print "Reading %s... " %(fdos), 
    sys.stdout.flush()
    doscar = open(fdos, 'r').readlines()
    sys.stdout.flush()

    self.natoms = int(doscar[0].strip().split(" ")[0])
    max_en, min_en, bins, fermi, x = [ float(e) for e in doscar[5].strip().split() ]
    self.bins = int(bins)
    _tdos = doscar[6:self.bins+6]
    ncol = len(_tdos[0].strip().split()) - 1
    if ncol == 4: self.spin = True
    en = np.zeros(self.bins)
    tdos = np.zeros((self.bins,ncol))
    for i, t in enumerate(_tdos):
      t = [ float(e) for e in t.strip().split() ]
      en[i] = t[0] - fermi
      for col in range(ncol):
        tdos[i][col] = t[col+1]/float(self.natoms)
    self.tdict = {}
    self.tdict['en'] = en
    self.tdict['tdos'] = tdos
    print "done."
    sys.stdout.flush()
    return
      
      
  def read_pdos(self, fdos='DOSCAR'):
    """
    Reads the projected density of states (DOS) data from the specified DOSCAR
    file, and stores the array of energies and DOS in a dictionary pdict. 
    In
    addition, calls the read_poscar() function to read element symbols.
    """
    if not os.path.exists(fdos):
      print "The file you specified does not exist."
      sys.exit(0)

    print "Reading %s... " %(os.path.abspath(fdos)),
    sys.stdout.flush()
    doscar = open(fdos, 'r').readlines()
    sys.stdout.flush()

    self.natoms = int(doscar[0].strip().split(" ")[0])
    max_en, min_en, bins, fermi, x = [ float(e) for e in doscar[5].strip().split() ]
    self.bins = int(bins)
    ntdos = self.bins+6
    _pdos = doscar[ntdos:]
    ncol = len(_pdos[1].strip().split()) - 1
    if ncol in [2, 8, 18]: self.spin = True
    en = np.zeros(self.bins)
    pdos = np.zeros((self.natoms,self.bins,ncol))
    for i, p in enumerate(_pdos[1:self.bins+1]):
      en[i] = float(p.strip().split()[0]) - fermi
    for i, p in enumerate(_pdos):
      if i%(self.bins+1) == 0: 
        nbin = 0
        continue
      site = i/(self.bins+1)
      p = [ float(e) for e in p.strip().split() ]
      for col in range(ncol):
        pdos[site][nbin][col] = p[col+1]
      nbin += 1
    self.pdict = {}
    self.pdict['en'] = en
    self.pdict['pdos'] = pdos
    print "done."
    sys.stdout.flush()
    return
    

  def read_params(self, fdos='DOSCAR', **kwargs):
    """
    Reads the input parameters from the specified input file. If none if
    specified, assumes that "params.in" resides in the same folder as the DOSCAR
    file.
    """
    dos_dir = os.path.dirname(os.path.abspath(fdos))
    fparams = os.path.join(dos_dir, 'params.in')
    if 'fparams' in kwargs:
      fparams = kwargs['fparams']
    print "Reading the input parameters from %s..." %(os.path.abspath(fparams)),
    params = open(fparams, 'r').readlines()
    self.projs = []
    [ self.projs.append(p.strip().split()) for p in params ]
    return


  def read_poscar(self, fdos='DOSCAR', **kwargs):
    """
    Reads the elements from the specified POSCAR/CONTCAR (fpos) file. If none is
    specified, assumes that the POSCAR file resides in the same directory as the
    DOSCAR file.
    """
    if 'fpos' in kwargs:
      fpos = kwargs['fpos']
    else:
      dos_dir = os.path.dirname(os.path.abspath(fdos))
      fpos = os.path.join(dos_dir, 'POSCAR')
    print "Reading the list of elements from %s..." %(fpos),
    pos = open(fpos, 'r').readlines()
    elems = pos[5].strip().split()
    if not all(isinstance(elem, str) for elem in elems):
      raise VaspFormatError("POSCAR/CONTCAR is not in VASP5 format.")
      sys.exit(0)
    natoms = [ int(n) for n in pos[6].strip().split() ]
    self._natoms = sum(natoms)
    self.elems = {}
    index = 1
    for elem, natom in zip(elems, natoms):
      for i in range(natom):
        self.elems[index] = elem
        index += 1
    print "done."
    return


  def plot_tdos(self, fdos='DOSCAR', **kwargs):
    """
    Plots the total density of states from the specified file into
    [fname][fformat] which defaults to tdos.jpg.
    """
    if self.tdict is None:
      self.read_tdos(fdos)
    en = self.tdict['en']
    dos = self.tdict['tdos']

    print "Now plotting... ",
    if self.spin:
      message = "Plotting spin-polarized DOS has not been implemented yet."
      raise ImplementationError(message)
      sys.exit(0)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(en, dos[:,0], color="#000000")
    ax.fill_between(en, 0, dos[:,0], color="#BBBBBB", alpha=0.5)
    ax.set_xlabel(r'$E - E_{\,\mathsf{Fermi}}$ (eV)', size=24)
    ax.set_ylabel(r'DOS (states/atom/eV)', size=24)
    if not 'plotrange' in kwargs:
      self._set_axes_limits(en, dos[:,0], ax)
    self._set_ticklabels(ax)

    fname = 'tdos'
    fformat = '.jpg'
    if 'fname' in kwargs: fname = kwargs['fname']
    if 'fformat' in kwargs: fformat = kwargs['fformat']
    dos_dir = os.path.dirname(os.path.abspath(fdos))
    ffig = os.path.join(dos_dir, fname+fformat)
    if 'plotrange' in kwargs: 
      plotrange = kwargs['plotrange']
      ffig = os.path.join(dos_dir, plotrange+'_'+fname+fformat)
    plt.savefig(ffig, bbox_inches='tight', dpi=300)
    print "done. tDOS plot in %s." %(ffig)
    return


  def plot_pdos(self, fdos='DOSCAR', total=True, **kwargs):
    """
    Plots the projected density of states from the specified file into
    [fname][fformat] which defaults to pdos.jpg.
    """
    if total:
      if self.tdict is None:
        self.read_tdos(fdos)
    if self.pdict is None:
      self.read_pdos(fdos)
    if self.elems is None:
      self.read_poscar(fdos)
    if self.projs is None:
      self.read_params(fdos)
    tdos = self.tdict['tdos']
    en = self.pdict['en']
    pdos = self.pdict['pdos']
    projs = self.projs

    print "Now plotting...",
    if self.spin:
      message = "Plotting spin-polarized DOS has not been implemented yet."
      raise ImplementationError(message)
      sys.exit(0)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    if total: 
      ax.plot(en, tdos[:,0], color='#bbbbbb', label="total")
      ax.fill_between(en, 0, tdos[:,0], color="#bbbbbb", alpha=0.40, label='total')
    for proj in projs:
      dos = self._calc_pdos_sum(pdos, proj)
      label = "%s-%s" %(proj[0], proj[3])
      color = proj[-1]
      ax.plot(en, dos, color=color, label=label)
      ax.fill_between(en, 0, dos, color=color, alpha=0.3)
      if total: continue
      if not 'plotrange' in kwargs:
        self._set_axes_limits(en, dos, ax)

    if not 'plotrange' in kwargs:
      if total: self._set_axes_limits(en, tdos[:,0], ax)
    ax.set_xlabel(r'$E - E_{\,\mathsf{Fermi}}$ (eV)', size=24)
    ax.set_ylabel(r'DOS (states/atom/eV)', size=24)
    if 'title' in kwargs:
      y1, y2 = ax.get_ylim()
      ax.text(-2, y2*0.9, r'%s' %(kwargs['title']), fontsize=28)
      #ax.set_title(r'%s' %(kwargs['title']), fontsize=28)
    self._set_ticklabels(ax)
    plt.legend(fontsize=20)

    fname = 'pdos'
    fformat = '.jpg'
    if 'fname' in kwargs: fname = kwargs['fname']
    if 'fformat' in kwargs: fformat = kwargs['fformat']
    dos_dir = os.path.dirname(os.path.abspath(fdos))
    ffig = os.path.join(dos_dir, fname+fformat)
    if 'plotrange' in kwargs: 
      plotrange = kwargs['plotrange']
      ffig = os.path.join(dos_dir, plotrange+'_'+fname+fformat)
    plt.savefig(ffig, bbox_inches='tight', dpi=300)
    print "done. pDOS plot in %s." %(ffig)
    

  def _set_axes_limits(self, en, dos, ax):
    if self.spin:
      message = "Plotting spin-polarized DOS has not been implemented yet."
      raise ImplementationError(message)
      sys.exit(0)
    

    xlim = [-6.0, 6.0]
    ax.set_xlim(xlim)
    ymax = 0
    for i in range(dos.size):
      if en[i] >= xlim[0] and en[i] <= xlim[1]:
        if dos[i] > ymax: ymax = dos[i]
    if ymax > 5.0:
      ymax = 5.0
    ylim = [ 0.0, 1.05*ymax ]
    ax.set_ylim(ylim)
    return


  def _set_ticklabels(self, ax):
    if self.spin:
      message = "Plotting spin-polarized DOS has not been implemented yet."
      raise ImplementationError(message)
      sys.exit(0)

    ax.xaxis.labelpad=10
    ax.yaxis.labelpad=10
    yloc = range(int(math.ceil(ax.get_ylim()[1])))
    ax.set_yticks(yloc)
    ax.set_yticklabels(yloc)
    for t in ax.xaxis.get_major_ticks():
      t.label.set_fontsize(20)
    for t in ax.yaxis.get_major_ticks():
      t.label.set_fontsize(20)
    return


  def _calc_pdos_sum(self, pdos, proj):
    elem = proj[0]
    sites = [ int(s) for s in proj[1].split('-') ]
    sites = range(sites[0]-1, sites[1])
    orbs = [ int(o) for o in proj[2].split('-') ]
    dos = np.zeros(self.bins)
    for site in sites:
      for orb in range(orbs[0]-1, orbs[1]):
        dos += pdos[site,:,orb]
    return dos/float(self.natoms)
