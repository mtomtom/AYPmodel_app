import numpy as np
from p2_functions.leafhoppers_class import Leafhoppers
from p2_functions.plants_class import Plants
# Random seed

class Simulation():
    def __init__(self):
        self.plants = None
        self.leafs = None
   
    def setUpPlants(self,fieldx,fieldy,infect,resolution,inf_plant1,inf_plant2,gamma1,gamma2,scaling,
    kGrow, mort_PI, grow_PU, timestep):
        self.plants = Plants(fieldx,fieldy,infect,resolution,inf_plant1,inf_plant2,gamma1,gamma2,scaling,
        kGrow, mort_PI, grow_PU, timestep)
        
    def setUpLeafhoppers(self, fieldx, fieldy,step,infect, uninf_leaf, inf_leaf1,inf_leaf2,resolution,
    gamma1,gamma2,scaling, rep, mort, ddmort, lh_direct, kl, n, timestep):
        self.leaf = Leafhoppers(fieldx, fieldy,step,infect,uninf_leaf, inf_leaf1, inf_leaf2,resolution,gamma1,gamma2,scaling,
        rep, mort, ddmort, lh_direct, kl, n, timestep)
        field_area = fieldx * fieldy
        init_uninf_leafhoppers = int(field_area * uninf_leaf*(1.0/scaling))
        init_inf_leafhoppers1 = int(field_area * inf_leaf1*(1.0/scaling))
        init_inf_leafhoppers2 = int(field_area * inf_leaf2*(1.0/scaling))

        self.leaf.inf_xcoord1 = np.random.rand(init_inf_leafhoppers1) * fieldx 
        self.leaf.inf_ycoord1 = np.random.rand(init_inf_leafhoppers1) * fieldy 
        self.leaf.inf_heading1  = np.random.rand(init_inf_leafhoppers1) * 2.0 * np.pi
        
        self.leaf.inf_xcoord2 = np.random.rand(init_inf_leafhoppers2) * fieldx
        self.leaf.inf_ycoord2 = np.random.rand(init_inf_leafhoppers2) * fieldy
        self.leaf.inf_heading2  = np.random.rand(init_inf_leafhoppers2) * 2.0 * np.pi

        self.leaf.uninf_xcoord = np.random.rand(init_uninf_leafhoppers) * fieldx
        self.leaf.uninf_ycoord = np.random.rand(init_uninf_leafhoppers) * fieldy
        self.leaf.uninf_heading  = np.random.rand(init_uninf_leafhoppers) * 2.0 * np.pi
               
    def moveLeafhoppers(self): 
        self.leaf.inf_xcoord1, self.leaf.inf_ycoord1, self.leaf.inf_heading1 = self.leaf.move(self.plants.matrix,self.leaf.inf_xcoord1, self.leaf.inf_ycoord1, self.leaf.inf_heading1)
        self.leaf.inf_xcoord2, self.leaf.inf_ycoord2, self.leaf.inf_heading2 = self.leaf.move(self.plants.matrix,self.leaf.inf_xcoord2, self.leaf.inf_ycoord2, self.leaf.inf_heading2)
        self.leaf.uninf_xcoord, self.leaf.uninf_ycoord, self.leaf.uninf_heading = self.leaf.move(self.plants.matrix,self.leaf.uninf_xcoord, self.leaf.uninf_ycoord, self.leaf.uninf_heading )
        
    def updateLeafhoppersPlants(self): 
        return self.plants.update_leafhoppers(self.leaf.inf_xcoord1,self.leaf.inf_ycoord1,self.leaf.inf_xcoord2,self.leaf.inf_ycoord2,self.leaf.uninf_xcoord,self.leaf.uninf_ycoord)   
    def infectLeafhoppers_density(self):
        return self.leaf.infectLeafhoppers_density(self.plants.matrix)
    def infectLeafhoppers_freq(self):
        return self.leaf.infectLeafhoppers_freq(self.plants.matrix)
    def infectLeafhoppers_bed(self):
        return self.leaf.infectLeafhoppers_bed(self.plants.matrix)
    def reproduction(self):
        return self.leaf.reproduceLeafhoppers(self.plants.matrix)
    def mortLeafhoppers(self):
        return self.leaf.mortLeafhoppers(self.plants.matrix)
    def updatePlants(self):
        return self.plants.update_plants()
    def addNewLeafs(self):
        return self.leaf.updateLeafhoppers()
    def shufflePlants(self):
        self.leaf.inf_xcoord1,self.leaf.inf_ycoord1,self.leaf.inf_xcoord2,self.leaf.inf_ycoord2,self.leaf.uninf_xcoord,self.leaf.uninf_ycoord = self.plants.shuffle_plants(self.leaf.inf_xcoord1,self.leaf.inf_ycoord1,self.leaf.inf_xcoord2,self.leaf.inf_ycoord2,self.leaf.uninf_xcoord,self.leaf.uninf_ycoord)
