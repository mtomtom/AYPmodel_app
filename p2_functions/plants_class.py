import numpy as np

class Plants():
    def __init__(self,fieldx,fieldy,infect,resolution,inf_plant1,inf_plant2,gamma1,gamma2,scaling,
    kGrow, mort_PI, grow_PU, timestep):
        self.x = float(fieldx) 
        self.y = float(fieldy)
        self.fieldx = float(fieldx)
        self.fieldy = float(fieldy)
        self.init_PHu = 1.0
        self.init_PHi1 = inf_plant1
        self.init_PHi2 = inf_plant2
        self.kGrow = kGrow
        self.infect = infect
        
        ## For the competition experiments, we scale virulence and transmission by a factor
        ## Transmission is only affected when applied from insects to plants!!!!
        
        # Scale the virulence by the severity factor
        self.mort_PI1 = mort_PI
        self.mort_PI2 = mort_PI
        self.grow_PU = grow_PU
        self.p_dict ={}
        # Different variables for different pathogen types
        self.p_dict["Hi1"]=0
        self.p_dict["Hi2"]=1
        self.p_dict["Hu"]=2
        self.p_dict["Li1"]=3
        self.p_dict["Li2"]=4
        self.p_dict["Lu"]=5
        self.p_dict["Crop"]=6

        self.timestep = timestep
        self.patches_x = int(self.x / resolution)
        self.patches_y = int(self.y / resolution)

        
        self.matrix = np.ones([int(self.patches_x),int(self.patches_y),6],dtype=np.float64)
        number_of_patches = self.patches_x * self.patches_y
        field_area = (self.x*self.y)
        self.patch_area = field_area / number_of_patches
        self.leafhopper_to_patch_x = self.patches_x / self.x
        self.leafhopper_to_patch_y = self.patches_y / self.y

        self.matrix[:,:,self.p_dict["Hu"]]=self.matrix[:,:,self.p_dict["Hu"]]*self.kGrow * self.init_PHu
        self.matrix[:,:,self.p_dict["Hi1"]]=self.matrix[:,:,self.p_dict["Hi1"]]*self.kGrow * self.init_PHi1
        self.matrix[:,:,self.p_dict["Hi2"]]=self.matrix[:,:,self.p_dict["Hi2"]]*self.kGrow * self.init_PHi2
        self.matrix[:,:,self.p_dict["Li1"]]*=0.0
        self.matrix[:,:,self.p_dict["Li2"]]*=0.0
        self.matrix[:,:,self.p_dict["Lu"]]*=0.0

        final_fields = np.zeros([self.patches_x,self.patches_y])
        self.matrix = np.dstack([self.matrix,final_fields])


    def update_leafhoppers(self,inf_xcoord1,inf_ycoord1,inf_xcoord2,inf_ycoord2,uninf_xcoord,uninf_ycoord):
        gridx = np.floor(inf_xcoord1 * self.leafhopper_to_patch_x).astype(int)
        gridy = np.floor(inf_ycoord1 * self.leafhopper_to_patch_y).astype(int)

        self.matrix[:,:,self.p_dict["Li1"]] = self.matrix[:,:,self.p_dict["Li1"]] * 0.0
        for i in range(len(gridx)):
            self.matrix[gridx[i],gridy[i],self.p_dict["Li1"]] +=(1.0 / self.patch_area)

        gridx = np.floor(inf_xcoord2 * self.leafhopper_to_patch_x).astype(int)
        gridy = np.floor(inf_ycoord2 * self.leafhopper_to_patch_y).astype(int)

        self.matrix[:,:,self.p_dict["Li2"]] = self.matrix[:,:,self.p_dict["Li2"]] * 0.0
        for i in range(len(gridx)):
            self.matrix[gridx[i],gridy[i],self.p_dict["Li2"]] +=(1.0 / self.patch_area)

        gridx = np.floor(uninf_xcoord * self.leafhopper_to_patch_x).astype(int)
        gridy = np.floor(uninf_ycoord * self.leafhopper_to_patch_y).astype(int)

        self.matrix[:,:,self.p_dict["Lu"]] = self.matrix[:,:,self.p_dict["Lu"]] * 0.0
        for i in range(len(gridx)):
            self.matrix[gridx[i],gridy[i],self.p_dict["Lu"]] +=(1.0 / self.patch_area)

    def update_plants(self):
        Hi1 = self.matrix[:, :,self.p_dict["Hi1"]]
        Hi2 = self.matrix[:, :,self.p_dict["Hi2"]]
        Hu = self.matrix[:, :,self.p_dict["Hu"]]
        Li1 = self.matrix[:, :,self.p_dict["Li1"]]
        Li2 = self.matrix[:, :,self.p_dict["Li2"]]
        
        def fHi(Hi,Hu,Li,HImort,alpha):
            return (alpha * Hu * Li) - (HImort * Hi)

        def fHu(Hi1,Hi2,Hu,Li1,Li2,HUgrow,alpha,kGrow):
            return (HUgrow * Hu * (1.0-((Hu + Hi1+Hi2) / kGrow))) - (alpha * Hu * (Li1+Li2))

        def rK3(Hi1,Hi2, Hu,Li1,Li2, fHi, fHu, timestep, HImort1,HImort2, alpha, HUgrow,kGrow):
            Hi1_1 = fHi(Hi1, Hu,Li1,HImort1,alpha)*timestep
            Hi2_1 = fHi(Hi2, Hu,Li2,HImort2,alpha)*timestep
            Hu1 = fHu(Hi1,Hi2, Hu,Li1,Li2,HUgrow,alpha,kGrow)*timestep

            Hi1k = Hi1 + Hi1_1*0.5
            Hi2k = Hi2 + Hi2_1*0.5
            Huk = Hu + Hu1*0.5

            Hi1_2 = fHi(Hi1k, Huk,Li1,HImort1,alpha)*timestep
            Hi2_2 = fHi(Hi2k, Huk,Li2,HImort2,alpha)*timestep
            Hu2 = fHu(Hi1k,Hi2k, Huk,Li1,Li2,HUgrow,alpha,kGrow)*timestep

            Hi1k = Hi1 + Hi1_2*0.5
            Hi2k = Hi2 + Hi2_2*0.5
            Hu2k = Hu + Hu2*0.5

            Hi1_3 = fHi(Hi1k, Hu2k,Li1,HImort1,alpha)*timestep
            Hi2_3 = fHi(Hi2k, Hu2k,Li2,HImort2,alpha)*timestep
            Hu3 = fHu(Hi1k,Hi2k, Hu2k,Li1,Li2,HUgrow,alpha,kGrow)*timestep

            Hi1k = Hi1 + Hi1_3
            Hi2k = Hi2 + Hi2_3
            Huk = Hu + Hu3

            Hi1_4 = fHi(Hi1k, Huk,Li1,HImort1,alpha)*timestep
            Hi2_4 = fHi(Hi2k, Huk,Li2,HImort2,alpha)*timestep
            Hu4 = fHu(Hi1k,Hi2k, Huk,Li1,Li2,HUgrow,alpha,kGrow)*timestep

            Hi1 = Hi1 + (Hi1_1 + 2*(Hi1_2 + Hi1_3) + Hi1_4)/6
            Hi2 = Hi2 + (Hi2_1 + 2*(Hi2_2 + Hi2_3) + Hi2_4)/6
            Hu = Hu + (Hu1 + 2*(Hu2 + Hu3) + Hu4)/6

            return Hi1,Hi2, Hu

        Hi1,Hi2, Hu= rK3(Hi1,Hi2, Hu, Li1,Li2, fHi, fHu, self.timestep, self.mort_PI1, self.mort_PI2,self.infect,self.grow_PU,self.kGrow)

        self.matrix[:, :,self.p_dict["Hi1"]]=Hi1
        self.matrix[:, :,self.p_dict["Hi2"]]=Hi2
        self.matrix[:, :,self.p_dict["Hu"]]=Hu

    def shuffle_plants(self,inf_xcoord1,inf_ycoord1,inf_xcoord2,inf_ycoord2,uninf_xcoord,uninf_ycoord):
        # Shuffle the plant matrix (For methods 1 and 2)
        # Shuffle each column separately, by the same index array
        # Create index array

        matrix,index = self.shuffle_columns(self.matrix,self.p_dict)
        # Update the infected and uninfected leafhopper matrices
        # First method: loop through every leafhopper and update its position
        reverse_index = np.empty(len(index))
        for i in range(len(index)):
            reverse_index[index[i]]=i
        # Reshape the reverse index array to match the plant matrix
        x_len = matrix.shape[0]
        y_len = matrix.shape[1]
        reverse_index.shape=(x_len,y_len)

        for i in range(len(inf_xcoord1)):
            grid_xcor = np.floor(inf_xcoord1[i] * self.leafhopper_to_patch_x).astype(int)
            grid_ycor = np.floor(inf_ycoord1[i] * self.leafhopper_to_patch_y).astype(int)

            new_coordinates_index = reverse_index[grid_xcor,grid_ycor].astype(int)
            new_coordinates = [np.floor(new_coordinates_index/y_len),np.floor(new_coordinates_index%y_len)]

            # Convert new patch coordinates into world coordinates
            inf_xcoord1[i] = (inf_xcoord1[i] - (grid_xcor/self.leafhopper_to_patch_x)) + (new_coordinates[0] / self.leafhopper_to_patch_x)
            inf_ycoord1[i] = (inf_ycoord1[i] - (grid_ycor/self.leafhopper_to_patch_y)) + (new_coordinates[1] / self.leafhopper_to_patch_y)

        for i in range(len(inf_xcoord2)):
            grid_xcor = np.floor(inf_xcoord2[i] * self.leafhopper_to_patch_x).astype(int)
            grid_ycor = np.floor(inf_ycoord2[i] * self.leafhopper_to_patch_y).astype(int)

            new_coordinates_index = reverse_index[grid_xcor,grid_ycor].astype(int)
            new_coordinates = [np.floor(new_coordinates_index/y_len),np.floor(new_coordinates_index%y_len)]

            # Convert new patch coordinates into world coordinates
            inf_xcoord2[i] = (inf_xcoord2[i] - (grid_xcor/self.leafhopper_to_patch_x)) + (new_coordinates[0] / self.leafhopper_to_patch_x)
            inf_ycoord2[i] = (inf_ycoord2[i] - (grid_ycor/self.leafhopper_to_patch_y)) + (new_coordinates[1] / self.leafhopper_to_patch_y)


        for i in range(len(uninf_xcoord)):
            grid_xcor = np.floor(uninf_xcoord[i] * self.leafhopper_to_patch_x).astype(int)
            grid_ycor = np.floor(uninf_ycoord[i] * self.leafhopper_to_patch_y).astype(int)

        #Method after talking to Stan
            new_coordinates_index = reverse_index[grid_xcor,grid_ycor].astype(int)
            new_coordinates = [np.floor(new_coordinates_index/y_len),np.floor(new_coordinates_index%y_len)]
            # Convert new patch coordinates into world coordinates
            uninf_xcoord[i] = (uninf_xcoord[i] - (grid_xcor/self.leafhopper_to_patch_x)) + (new_coordinates[0] / self.leafhopper_to_patch_x)
            uninf_ycoord[i] = (uninf_ycoord[i] - (grid_ycor/self.leafhopper_to_patch_y)) + (new_coordinates[1] / self.leafhopper_to_patch_y)

            self.matrix = matrix

        return (inf_xcoord1,inf_ycoord1,inf_xcoord2,inf_ycoord2,uninf_xcoord,uninf_ycoord)


    def shuffle_columns(self,matrix,p_dict):
        xdim = matrix.shape[0]
        ydim = matrix.shape[1]
        index = np.arange(xdim * ydim)
        np.random.shuffle(index)
        # Shuffle each column of the plants matrix, according to the index

        temp = matrix[:,:,p_dict["Hi1"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Hi1"]] = temp

        temp = matrix[:,:,p_dict["Hi2"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Hi2"]] = temp

        temp = matrix[:,:,p_dict["Hu"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Hu"]] = temp

        temp = matrix[:,:,p_dict["Li1"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Li1"]] = temp

        temp = matrix[:,:,p_dict["Li2"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Li2"]] = temp

        temp = matrix[:,:,p_dict["Lu"]]
        temp.shape = (xdim * ydim)
        temp = temp[index]
        temp.shape=(xdim,ydim)
        matrix[:,:,p_dict["Lu"]] = temp

        return matrix,index
