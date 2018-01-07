#WillkÃ¼rliche Festlegung: Gausspaket soll Breite der Barriere/5 breit sein


import numpy as np

class Potential():
    """
    Ceates an potential that is suitable for the following calculations, by
        prefixing an relevant amount of zeros in front of the given potential.
        It also appends a high amount behind the given Potential to keep the wave monitored.
        
        Creates a position grid to be used in the Fourier- as well as SplitStep-class.
        
        Calculates the position to be used for the gaussian wave symmetry point.
        
        NOTICE: barrier.size should be 20 or greater to provide an satisfying resoluion for the gauss-package
    
    Parameters
    ----------
    barrier : Array
        Potentialvalues on discret samplepoint-grid in position-space
    dx : Float
        distance between two neighbouring samplepoints in position space grid
    
    Attributes
    ----------
    position_grid : Array
        grid of considered positions/ sample points in position-space
    pos_grid_width : float
        total width of the potential
    potential: Array
        useable potential to be used for futher calculaions
    gauss_index_width : integer
        Number of samplepoints corresponding to the gaussian package width
    gauss_width : float
        width to be used for the gaussian package in position space
    gauss_symmetry_index : integer
        self.position_grid[self.gauss_symmetry_index] is the point in position grid to be used for the gaussian symmetrypoint
    """
    
    def __init__(self, barrier, dx):
        
        self.barrier = barrier
        self.dx = dx
           
        def create_gauss_index_width(self):
            barrier_index_width = self.barrier.size
            return int(barrier_index_width/5)
        
        self.gauss_index_width = create_gauss_index_width(self)
        self.gauss_width = self.gauss_index_width*self.dx
        
        def create_potential(self):
            praefix = np.zeros(int(6*self.gauss_index_width+1))
            suffix = np.zeros(500*self.gauss_index_width)
            
            pot = np.append(praefix,self.barrier)
            return np.append(pot, suffix)
            
        self.potential = create_potential(self)
        
        def create_position_grid(self):
            potential_index_width = self.potential.size
            return np.linspace(0,(potential_index_width)*self.dx, potential_index_width)
        
        self.position_grid = create_position_grid(self)
        
        def create_posistion_grid_width(self):
            return self.position_grid[-1]-self.position_grid[0]
        
        self.pos_grid_width = create_posistion_grid_width(self)
        
        def create_gauss_symmetry_index(self):
            return 3*self.gauss_index_width
        
        self.gauss_symmetry_index = create_gauss_symmetry_index(self)
        self.gauss_symmerey_point = self.position_grid[self.gauss_symmetry_index]
