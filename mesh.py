import numpy as np

# ---------------------
# building mesh arrays 
# ---------------------

class Mesh():
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'                
        
        # -----
        # Grid
        # -----
        with open(directory+'radial_grid', 'r') as f:
            self.ndomains = int(f.readline())
            self.nzd = np.zeros(self.ndomains, dtype=int)
            for i in range(self.ndomains):
                values = f.readline().split()
                self.nzd[i] = int(values[1])
            self.r = self.read_scalar(f, self.nr + 1)

        with open(directory+'lat_grid', 'r') as f:
            self.theta = self.read_scalar(f, self.nth)

        with open(directory+'2D_grid', 'r') as f:
            f.readline()
            self.x = self.read_array(f, self.nr+1, self.nth)
            f.readline()
            self.y = self.read_array(f, self.nr+1, self.nth)
