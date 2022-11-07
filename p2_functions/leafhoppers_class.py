import numpy as np
from scipy import signal
from scipy.ndimage.filters import gaussian_filter
from bisect import bisect

class Leafhoppers():
    def __init__(self,fieldx, fieldy, step,infect, uninf_leaf, inf_leaf1, inf_leaf2,resolution,gamma1,gamma2,
    scaling, rep, mort, ddmort, lh_direct, kl, n, timestep):
        # Parameters
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.fieldx = fieldx
        self.fieldy = fieldy
        self.step = step
        self.uninf_rep_Hi = rep
        self.uninf_rep_Hu = rep
        self.inf_rep_Hi = rep
        self.inf_rep_Hu = rep
        self.infect = infect
        self.mort = mort
        self.ddmort =ddmort
        self.lh_direct = lh_direct
        self.inf_init_conc1 = inf_leaf1
        self.inf_init_conc2 = inf_leaf2
        self.uninf_init_conc = uninf_leaf
        self.kl = kl
        self.n = n
        self.timestep = timestep
        self.patches_x = int(self.fieldx / resolution)
        self.patches_y = int(self.fieldy / resolution)
        self.lhop_to_patch_x = self.patches_x / self.fieldx
        self.lhop_to_patch_y = self.patches_y / self.fieldy
        self.field_area = fieldx * fieldy
        self.p_dict ={}
        # Different variables for different pathogen types
        self.p_dict["Hi1"]=0
        self.p_dict["Hi2"]=1
        self.p_dict["Hu"]=2
        self.p_dict["Li1"]=3
        self.p_dict["Li2"]=4
        self.p_dict["Lu"]=5
        
        # Variables
        self.inf_xcoord1 = None
        self.inf_ycoord1 = None
        self.inf_xcoord2 = None
        self.inf_ycoord2 = None
        self.inf_heading1 = None
        self.inf_heading2 = None
        self.uninf_xcoord = None
        self.uninf_ycoord = None
        self.uninf_heading = None
        
                
    def move(self,plants,xcoord,ycoord,heading):
        # Infected leafhoppers       
        inf_leaf_no = len(xcoord)       
        # Creates an array of angles (in radians) - only rotate once a day
        elements_to_change = inf_leaf_no * self.timestep * self.lh_direct # prob of leafhoppers moving - change direction 10x a day
    ###############################################################################################
    # If time step * number of leafhoppers falls below 1, then this method will fail
    # Add in fix to stop this happening
    #################################################################################################
        if elements_to_change < 1:
            elements_to_change = 1
        

        leafhopper_step = self.step * self.timestep
    
        mask = np.random.randint(0,elements_to_change,size=heading.shape).astype(np.bool)
        new_heading = np.random.rand(inf_leaf_no) * 2.0 * np.pi
        # code for a biased random walk - moving north 
        heading[mask] = new_heading[mask]
        heading = np.random.rand(inf_leaf_no)*2.0*np.pi
        xcoord += (leafhopper_step * np.cos(heading))
        ycoord += (leafhopper_step * np.sin(heading))
        xcoord[np.where(xcoord >= self.fieldx)] -= self.fieldx
        xcoord[np.where(xcoord < 0.0)] += self.fieldx
        ycoord[np.where(ycoord >= self.fieldy)] -= self.fieldy
        ycoord[np.where(ycoord < 0.0)] += self.fieldy
        return xcoord, ycoord, heading

        
    def infectLeafhoppers_freq(self,plants):
        ## Frequency dependent transmission version of the function
        # Function to choose between no infection, infection1, infection2

        def weighted_choice(choices):
            values, weights = zip(*choices)
            total = 0
            cum_weights = []
            for w in weights:
                total += w
                cum_weights.append(total)
                
            x = np.random.rand() * total
            i = bisect(cum_weights, x)
            return values[i]
        
        def select_infect(inf1,inf2,uninf):
            weights = np.array([inf1,inf2,uninf])
            rw = 1.0
            while rw > 0:
                index = np.random.randint(0,3,1)
                rw = np.random.rand()
                rw = rw - weights[index]
            
            return index
        
        #########################

        # Create infection matrix
        gridx = np.floor(self.uninf_xcoord * self.lhop_to_patch_x).astype(int)
        gridy = np.floor(self.uninf_ycoord * self.lhop_to_patch_y).astype(int)
        uninf_leaf_no = len(self.uninf_xcoord)
        inf_leaf_no = len(self.inf_heading1)
   
        infect_dict = {}
        infect1 = self.infect * self.gamma1
        infect2 = self.infect * self.gamma2

        # Create a value for each leafhopper, based on the plants in its grid cell

        total_plants = plants[gridx,gridy,self.p_dict["Hi1"]] + plants[gridx,gridy,self.p_dict["Hi2"]] + plants[gridx,gridy,self.p_dict["Hu"]]
   
        infection_array1 = infect1 * self.timestep * (plants[gridx,gridy,self.p_dict["Hi1"]] / total_plants)
        infection_array2 = infect2 * self.timestep * (plants[gridx,gridy,self.p_dict["Hi2"]] / total_plants)
 
        result = np.zeros(uninf_leaf_no)
        for i in range(len(result)):
            inf1 = infection_array1[i]
            inf2 = infection_array2[i]
            uninf = 1.0-(inf1+inf2)
            result[i]=weighted_choice([(0,uninf), (1,inf1), (2,inf2)])
            #result[i]=select_infect(uninf,inf1,inf2)
            
        infection_mask1= result == 1
        infection_mask2= result == 2
        
        #random_array1 = np.random.rand(uninf_leaf_no)
        #random_array2 = np.random.rand(uninf_leaf_no)
        
        #infection_mask1 = random_array1 < infection_array1
        #infection_mask2 = random_array2 < infection_array2

        self.new_infected_x1 = self.uninf_xcoord[infection_mask1]
        self.new_infected_y1 = self.uninf_ycoord[infection_mask1]
        self.new_infected_heading1 = self.uninf_heading[infection_mask1]
        self.new_infected_x2 = self.uninf_xcoord[infection_mask2]
        self.new_infected_y2 = self.uninf_ycoord[infection_mask2]
        self.new_infected_heading2 = self.uninf_heading[infection_mask2]
        
        total_infection = [any(tup) for tup in zip(infection_mask1, infection_mask2)]
        total_infection = np.array(total_infection).astype(bool)
        
        self.uninf_xcoord = self.uninf_xcoord[~total_infection]
        self.uninf_ycoord = self.uninf_ycoord[~total_infection]
        self.uninf_heading = self.uninf_heading[~total_infection]
        
        self.inf_xcoord1 = np.hstack([self.inf_xcoord1,self.new_infected_x1])
        self.inf_ycoord1 = np.hstack([self.inf_ycoord1,self.new_infected_y1])
        self.inf_heading1 = np.hstack([self.inf_heading1,self.new_infected_heading1])
        self.inf_xcoord2 = np.hstack([self.inf_xcoord2,self.new_infected_x2])
        self.inf_ycoord2 = np.hstack([self.inf_ycoord2,self.new_infected_y2])
        self.inf_heading2 = np.hstack([self.inf_heading2,self.new_infected_heading2])

    def infectLeafhoppers_density(self,plants):
        # Function to choose between no infection, infection1, infection2

        def weighted_choice(choices):
            values, weights = zip(*choices)
            total = 0
            cum_weights = []
            for w in weights:
                total += w
                cum_weights.append(total)
                
            x = np.random.rand() * total
            i = bisect(cum_weights, x)
            return values[i]
        
        def select_infect(inf1,inf2,uninf):
            weights = np.array([inf1,inf2,uninf])
            rw = 1.0
            while rw > 0:
                index = np.random.randint(0,3,1)
                rw = np.random.rand()
                rw = rw - weights[index]
            
            return index

        # Create infection matrix
        gridx = np.floor(self.uninf_xcoord * self.lhop_to_patch_x).astype(int)
        gridy = np.floor(self.uninf_ycoord * self.lhop_to_patch_y).astype(int)
        uninf_leaf_no = len(self.uninf_xcoord)
        
        infect_dict = {}
        infect1 = self.infect * self.gamma1
        infect2 = self.infect * self.gamma2

        # Create a value for each leafhopper, based on the plants in its grid cell
   
        infection_array1 = infect1 * self.timestep * plants[gridx,gridy,self.p_dict["Hi1"]]
        infection_array2 = infect2 * self.timestep * plants[gridx,gridy,self.p_dict["Hi2"]]
        
        result = np.zeros(uninf_leaf_no)
        for i in range(len(result)):
            inf1 = infection_array1[i]
            inf2 = infection_array2[i]
            uninf = 1.0-(inf1+inf2)
            result[i]=weighted_choice([(0,uninf), (1,inf1), (2,inf2)])
            #result[i]=select_infect(uninf,inf1,inf2)
            
        infection_mask1= result == 1
        infection_mask2= result == 2
        
        #random_array1 = np.random.rand(uninf_leaf_no)
        #random_array2 = np.random.rand(uninf_leaf_no)
        
        #infection_mask1 = random_array1 < infection_array1
        #infection_mask2 = random_array2 < infection_array2

        self.new_infected_x1 = self.uninf_xcoord[infection_mask1]
        self.new_infected_y1 = self.uninf_ycoord[infection_mask1]
        self.new_infected_heading1 = self.uninf_heading[infection_mask1]
        self.new_infected_x2 = self.uninf_xcoord[infection_mask2]
        self.new_infected_y2 = self.uninf_ycoord[infection_mask2]
        self.new_infected_heading2 = self.uninf_heading[infection_mask2]
        
        total_infection = [any(tup) for tup in zip(infection_mask1, infection_mask2)]
        total_infection = np.array(total_infection).astype(bool)
        
        self.uninf_xcoord = self.uninf_xcoord[~total_infection]
        self.uninf_ycoord = self.uninf_ycoord[~total_infection]
        self.uninf_heading = self.uninf_heading[~total_infection]
        
        self.inf_xcoord1 = np.hstack([self.inf_xcoord1,self.new_infected_x1])
        self.inf_ycoord1 = np.hstack([self.inf_ycoord1,self.new_infected_y1])
        self.inf_heading1 = np.hstack([self.inf_heading1,self.new_infected_heading1])
        self.inf_xcoord2 = np.hstack([self.inf_xcoord2,self.new_infected_x2])
        self.inf_ycoord2 = np.hstack([self.inf_ycoord2,self.new_infected_y2])
        self.inf_heading2 = np.hstack([self.inf_heading2,self.new_infected_heading2])
 
    def reproduceLeafhoppers(self,plants):
        gridx = np.floor(self.uninf_xcoord * self.lhop_to_patch_x).astype(int)
        gridy = np.floor(self.uninf_ycoord * self.lhop_to_patch_y).astype(int)
        uninf_leaf_no = len(self.uninf_xcoord)
        # Create a value for each leafhopper, based on the plants in its grid cell
        reproduce_array = self.uninf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi1"]])
        reproduce_array = reproduce_array + (self.uninf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi2"]]))
        reproduce_array = reproduce_array + (self.uninf_rep_Hu * self.timestep *
                                         (plants[gridx,gridy,self.p_dict["Hu"]]))
        # Uninfected leafhopper reproduction
        random_array = np.random.rand(uninf_leaf_no)
        reproduce_array_mask = random_array < reproduce_array
        self.new_uninf_xcoord = self.uninf_xcoord[reproduce_array_mask]
        self.new_uninf_ycoord = self.uninf_ycoord[reproduce_array_mask]
        new_heading = self.uninf_heading[reproduce_array_mask]
        self.new_uninf_heading = np.random.rand(len(new_heading)) * 2.0 * np.pi
        
        # Infected1 leafhopper reproduction
        gridx = np.floor(self.inf_xcoord1 * self.lhop_to_patch_x).astype(int)
        gridy = np.floor(self.inf_ycoord1 * self.lhop_to_patch_y).astype(int)
        inf_leaf_no = len(self.inf_xcoord1)
        # Create a value for each leafhopper, based on the plants in its grid cell
        reproduce_array = self.inf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi1"]])
        reproduce_array = reproduce_array + (self.inf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi2"]]))
        reproduce_array = reproduce_array + (self.inf_rep_Hu * self.timestep *
                                         (plants[gridx,gridy,self.p_dict["Hu"]]))

        random_array = np.random.rand(inf_leaf_no)
        reproduce_array_mask = random_array < reproduce_array
        self.new_inf_xcoord1 = self.inf_xcoord1[reproduce_array_mask]
        self.new_inf_ycoord1 = self.inf_ycoord1[reproduce_array_mask]
        new_heading1 = self.inf_heading1[reproduce_array_mask]
        self.new_inf_heading1 = np.random.rand(len(new_heading1)) * 2.0 * np.pi
        
        # Infected2 leafhopper reproduction
        gridx = np.floor(self.inf_xcoord2 * self.lhop_to_patch_x).astype(int)
        gridy = np.floor(self.inf_ycoord2 * self.lhop_to_patch_y).astype(int)
        inf_leaf_no = len(self.inf_xcoord2)
        # Create a value for each leafhopper, based on the plants in its grid cell
        reproduce_array = self.inf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi1"]])
        reproduce_array = reproduce_array + (self.inf_rep_Hi * self.timestep * (plants[gridx,gridy,self.p_dict["Hi2"]]))
        reproduce_array = reproduce_array + (self.inf_rep_Hu * self.timestep *
                                         (plants[gridx,gridy,self.p_dict["Hu"]]))
        random_array = np.random.rand(inf_leaf_no)
        reproduce_array_mask = random_array < reproduce_array
        self.new_inf_xcoord2 = self.inf_xcoord2[reproduce_array_mask]
        self.new_inf_ycoord2 = self.inf_ycoord2[reproduce_array_mask]
        new_heading2 = self.inf_heading2[reproduce_array_mask]
        self.new_inf_heading2 = np.random.rand(len(new_heading2)) * 2.0 * np.pi
        
    def mortLeafhoppers(self,plants):
        uninf_leaf_no = len(self.uninf_xcoord)
        inf_leaf_no1 = len(self.inf_xcoord1)
        inf_leaf_no2 = len(self.inf_xcoord2)
        # Linear mortality rate
        death_array_uninf = np.ones(uninf_leaf_no) * self.mort * self.timestep
        death_array_inf1 = np.ones(inf_leaf_no1) * self.mort * self.timestep
        death_array_inf2 = np.ones(inf_leaf_no2) * self.mort * self.timestep
        
        rand_array_uninf = np.random.rand(uninf_leaf_no)
        rand_array_inf1 = np.random.rand(inf_leaf_no1)
        rand_array_inf2 = np.random.rand(inf_leaf_no2)
        
        death_list_mask_uninf = rand_array_uninf < death_array_uninf
        death_list_mask_inf1 = rand_array_inf1 < death_array_inf1
        death_list_mask_inf2 = rand_array_inf2 < death_array_inf2
        # Remove dead leafhoppers
        self.uninf_xcoord = self.uninf_xcoord[~death_list_mask_uninf]
        self.uninf_ycoord = self.uninf_ycoord[~death_list_mask_uninf]
        self.uninf_heading = self.uninf_heading[~death_list_mask_uninf]
        
        self.inf_xcoord1 = self.inf_xcoord1[~death_list_mask_inf1]
        self.inf_ycoord1 = self.inf_ycoord1[~death_list_mask_inf1]
        self.inf_heading1 = self.inf_heading1[~death_list_mask_inf1]
        
        self.inf_xcoord2 = self.inf_xcoord2[~death_list_mask_inf2]
        self.inf_ycoord2 = self.inf_ycoord2[~death_list_mask_inf2]
        self.inf_heading2 = self.inf_heading2[~death_list_mask_inf2]

        uninf_leaf_no = len(self.uninf_xcoord)
        inf_leaf_no1 = len(self.inf_xcoord1)
        inf_leaf_no2 = len(self.inf_xcoord2)
        
        total_l_conc = (uninf_leaf_no + inf_leaf_no1 + inf_leaf_no2) / self.field_area
        death_array_uninf = np.ones(uninf_leaf_no) * self.ddmort * self.timestep * ((total_l_conc**self.n)/((total_l_conc**self.n)+(self.kl**self.n)))
        death_array_inf1 = np.ones(inf_leaf_no1) * self.ddmort * self.timestep * ((total_l_conc**self.n)/((total_l_conc**self.n)+(self.kl**self.n)))
        death_array_inf2 = np.ones(inf_leaf_no2) * self.ddmort * self.timestep * ((total_l_conc**self.n)/((total_l_conc**self.n)+(self.kl**self.n)))
    
        # create grid array
        grid_xcor_uninf = np.floor(self.uninf_xcoord* self.lhop_to_patch_x).astype(int)
        grid_ycor_uninf = np.floor(self.uninf_ycoord* self.lhop_to_patch_y).astype(int)
        grid_xcor_inf1 = np.floor(self.inf_xcoord1* self.lhop_to_patch_x).astype(int)
        grid_ycor_inf1 = np.floor(self.inf_ycoord1* self.lhop_to_patch_y).astype(int)
        grid_xcor_inf2 = np.floor(self.inf_xcoord2* self.lhop_to_patch_x).astype(int)
        grid_ycor_inf2 = np.floor(self.inf_ycoord2* self.lhop_to_patch_y).astype(int)
    
        plants_average = plants[:,:,self.p_dict["Li1"]] + plants[:,:,self.p_dict["Li2"]]+ plants[:,:,self.p_dict["Lu"]]
        
        # !!!!!!!!!! MAKE SURE YOU CHANGE THE BLUR RADIUS WHEN YOU CHANGE THE PATCH SIZE !!!!!!!!!!!
        # default = 10.0
        sigma_scaled = 10.0 * self.lhop_to_patch_x
        conv_out = gaussian_filter(plants_average, sigma_scaled,mode="wrap")

        death_array_uninf =  np.ones(uninf_leaf_no) * self.ddmort * self.timestep *(
                                        ((conv_out[grid_xcor_uninf,grid_ycor_uninf])**self.n)/
                                        (((conv_out[grid_xcor_uninf,grid_ycor_uninf])**self.n)+self.kl**self.n))
        
        death_array_inf1 =  np.ones(inf_leaf_no1) * self.ddmort * self.timestep *(
                                        ((conv_out[grid_xcor_inf1,grid_ycor_inf1])**self.n)/
                                        (((conv_out[grid_xcor_inf1,grid_ycor_inf1])**self.n)+self.kl**self.n))
        
        death_array_inf2 =  np.ones(inf_leaf_no2) * self.ddmort * self.timestep *(
                                        ((conv_out[grid_xcor_inf2,grid_ycor_inf2])**self.n)/
                                        (((conv_out[grid_xcor_inf2,grid_ycor_inf2])**self.n)+self.kl**self.n))
        
        rand_array_uninf = np.random.rand(uninf_leaf_no)
        rand_array_inf1 = np.random.rand(inf_leaf_no1)
        rand_array_inf2 = np.random.rand(inf_leaf_no2)
        death_list_mask_uninf = rand_array_uninf < death_array_uninf
        death_list_mask_inf1 = rand_array_inf1 < death_array_inf1
        death_list_mask_inf2 = rand_array_inf2 < death_array_inf2
        # Remove dead leafhoppers
        self.uninf_xcoord = self.uninf_xcoord[~death_list_mask_uninf]
        self.uninf_ycoord = self.uninf_ycoord[~death_list_mask_uninf]
        self.uninf_heading = self.uninf_heading[~death_list_mask_uninf]
        
        self.inf_xcoord1 = self.inf_xcoord1[~death_list_mask_inf1]
        self.inf_ycoord1 = self.inf_ycoord1[~death_list_mask_inf1]
        self.inf_heading1 = self.inf_heading1[~death_list_mask_inf1]
        
        self.inf_xcoord2 = self.inf_xcoord2[~death_list_mask_inf2]
        self.inf_ycoord2 = self.inf_ycoord2[~death_list_mask_inf2]
        self.inf_heading2 = self.inf_heading2[~death_list_mask_inf2]
        
    def updateLeafhoppers(self):
        # Add the new leafhoppers to the array
        if self.new_inf_xcoord1.any():
            self.inf_xcoord1 = np.hstack([self.inf_xcoord1,self.new_inf_xcoord1])
            self.inf_ycoord1 = np.hstack([self.inf_ycoord1,self.new_inf_ycoord1])
            self.inf_heading1 = np.hstack([self.inf_heading1,self.new_inf_heading1])
            self.new_inf_xcoord1 = 0
            self.new_inf_ycoord1 = 0
            self.new_inf_heading1 = 0
            
        if self.new_inf_xcoord2.any():
            self.inf_xcoord2 = np.hstack([self.inf_xcoord2,self.new_inf_xcoord2])
            self.inf_ycoord2 = np.hstack([self.inf_ycoord2,self.new_inf_ycoord2])
            self.inf_heading2 = np.hstack([self.inf_heading2,self.new_inf_heading2])
            self.new_inf_xcoord2 = 0
            self.new_inf_ycoord2 = 0
            self.new_inf_heading2 = 0
    
        if self.new_uninf_xcoord.any():
            self.uninf_xcoord = np.hstack([self.uninf_xcoord,self.new_uninf_xcoord])
            self.uninf_ycoord = np.hstack([self.uninf_ycoord,self.new_uninf_ycoord])
            self.uninf_heading = np.hstack([self.uninf_heading,self.new_uninf_heading])
            self.new_uninf_xcoord = 0
            self.new_uninf_ycoord = 0
            self.new_uninf_heading = 0

