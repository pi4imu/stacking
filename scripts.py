# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r, centroid=True, delete_bright_regions=True, delete_superfluous=True, draw=True, draw_additional=False, redshifted_back=True):

    # there are several cases of SAME ihal for DIFFERENT cluster numbers
    # this is the reason for using cluster number as a counter
    
    current_cluster = clusters.loc[current_cluster_number]
    
    RA_c = current_cluster["x_pix"]*30-5
    DEC_c = current_cluster["y_pix"]*30-5
    R_vir = current_cluster["Rrel"]*30
    R_500 = current_cluster["R500"]*0.704  # kpc
    ztrue = current_cluster["z_true"]
    
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    
    R_500_rescaled = R_500/D_A.value*180/np.pi # degrees
    
    D_A = 343000
    
    R_500_fid = 1000/343000*180/np.pi   # degrees
    
    snap_id_str = binned_clusters[current_cluster_number][1]   # id of photon list
        
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
      
    #t1 = Table.read("../data/eROSITA_30.0x30.0/Phox/AGNphlist_"+snap_id_str+".fits", hdu=2)
    
    #t = vstack([t, t1])
    
    SLICE = t.to_pandas()        # for photons extraction
    SLICE1 = t.to_pandas()       # for drawing
    SLICE2 = t.to_pandas()       # for center searching if it is not from table
    SLICE3 = t.to_pandas()       # for rescaling SLICE2
    
    R = r * R_500_fid
    
    print(R, R*60)
        
    AREA = np.pi*R**2*3600   # min2
        
    if not centroid:
    
        # taking photons from circle centered at RA_c, DEC_c
    
        SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        cntr = (RA_c, DEC_c) # for drawing
    
    else:
    
        #setting area and resolution for searching for center
    
        ang_res = 4
        halfsidelength = 5                    # in R500
        half_size = halfsidelength*R   # in degrees
        
        if (current_cluster_number != 13334) and (current_cluster_number != 18589):
            hs4s = half_size/3/(halfsidelength/3)
        else:
            hs4s = half_size/1/(halfsidelength/3)
            
        # making 2D histogram with side length 2*half_size with center (RA_c, DEC_c) without drawing
                
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < hs4s) & (np.abs(SLICE2["DEC"]-DEC_c) < hs4s), True, False)
        whattodraw = SLICE2[SLICE2['what'] == True]
        whattodraw = whattodraw.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], bins=int(2*hs4s*3600/ang_res))
                       
        # centroid position
                
        psum = sum(nmhg.flatten())
        c_x, c_y = 0, 0
        
        for i in range(0, len(nmhg)):
            for j in range(0, len(nmhg)):
                c_x = c_x + i*nmhg[i,j]/psum
                c_y = c_y + j*nmhg[i,j]/psum
        
        # position of centroid in units of pixels relative to the upper left border    
        c = [int(c_x), len(nmhg)-int(c_y)]
        
        if draw_additional:
            plt.imshow(np.rot90(nmhg),
                       norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                       origin='upper')
            plt.scatter(c[0], c[1], color='red')
            plt.scatter(int(len(nmhg)/2), int(len(nmhg)/2), color='magenta')
            plt.gca().set_aspect('equal', 'box')
            plt.xlim(c[0]-int(hs4s*3600/ang_res), c[0]+int(hs4s*3600/ang_res))
            plt.ylim(c[1]+int(hs4s*3600/ang_res), c[1]-int(hs4s*3600/ang_res))
            plt.gca().add_patch(plt.Circle(c, int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))
            plt.show()
                
        c_x_1 = RA_c - hs4s + c_x*ang_res/3600 
        c_y_1 = DEC_c - hs4s + c_y*ang_res/3600
        cntr = (c_x_1, c_y_1) # position in degrees
        
        # taking photons from circle centered at centroid
        
        SLICE["check"]=np.where((SLICE["RA"]-c_x_1)**2+(SLICE["DEC"]-c_y_1)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        # deleting less massive haloes
        
        if delete_superfluous:
            
            VICINITY = np.where( ((clusters_all["x_pix"]*30-5 - c_x_1)**2 + (clusters_all["y_pix"]*30-5 - c_y_1)**2 < half_size**2) & ( np.abs(clusters_all["z_true"] - ztrue) < 0.017141 ) )
            
            #print(VICINITY)
            
            vclu = clusters_all.loc[VICINITY]
            
            for clcl in VICINITY:
            
                vicinity_current = clusters_all.loc[clcl]
               
                vicenter = list(zip(vicinity_current["x_pix"].values*30-5, vicinity_current["y_pix"].values*30-5))
        
                #print(vicenter)
        
        # deleting bright regions           
        
        if delete_bright_regions:
        
            number_of_unfiltered_photons = len(dddfff)
        
            # recalculate nmhg relative to centroid
            
            SLICE3["what"] = np.where( (np.abs(SLICE3["RA"]-c_x_1) < half_size) & (np.abs(SLICE3["DEC"]-c_y_1) < half_size),
                                      True, False)
            whattodraw = SLICE3[SLICE3['what'] == True]
            whattodraw = whattodraw.drop("what", axis=1)
            nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], 
                                        bins=int(2*half_size*3600/ang_res),
                                        range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                        (cntr[1]-half_size, cntr[1]+half_size)]))            
            
            #fltr = Tophat2DKernel(2)
            fltr = Gaussian2DKernel(1)
            nmhg = convolve_fft(nmhg, fltr) 

            # some magic ...
            
            #shift = [int((RA_c-cntr[0]+half_size)*3600/ang_res), int((DEC_c-cntr[1]+half_size)*3600/ang_res)]
            shift = [int((half_size)*3600/ang_res), int((half_size)*3600/ang_res)]
            c1 = shift 
                       
            nmhg1 = kruzhok(int(R*3600/ang_res), c1, nmhg, int(R*3600/ang_res)+1)[0]
  
            if draw_additional:
                
                plt.show()
            
                # initial plot in terms of pixels and only circle of R500
        
                plt.figure(figsize=(12,5))
                plt.subplot(121)
                plt.title("nmhg")
                plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.scatter(c1[0], c1[1], color='orangered')
                plt.gca().add_patch(plt.Circle(c1, int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.subplot(122)
                plt.title("nmhg1, R = "+str(int(R*3600/ang_res))+" pixels")
                plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.gca().add_patch(plt.Circle((int(R*3600/ang_res), int(R*3600/ang_res)), 
                                    int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))  
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.show()
            
            # searching for cutoff value in one pixel            
            
            threshold = 0.01
            
            beeens = np.geomspace(1, np.max(nmhg1.flatten()), 50)
            amount_in_bin, bin_borders = np.histogram(nmhg1.flatten(), bins = beeens)
            bin_centers = (bin_borders[:-1]+bin_borders[1:])/2
            cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))]
            number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
            index_cutoff = list(bin_centers).index(number_cutoff)
            threshold_new = sum(amount_in_bin[:index_cutoff]) / cumulative[-1]
           
            if draw_additional:
            
                plt.figure(figsize=(13,5))
                plt.subplot(121)
                plt.hist(nmhg1.flatten(), bins = beeens)
                plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
                plt.ylabel("Amount of such bins")
                plt.yscale("log")
                plt.xscale("log")
                plt.title("Flattened histogram for upper right image")
                plt.subplot(122)
                plt.scatter(bin_centers, cumulative)
                plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
                plt.ylabel("Cumulative distribution")
                plt.yscale("log")
                plt.xscale("log")         
                plt.title(f"$x$-values are added up until their sum\nis right below {threshold*100:.0f} % cutoff")    
                plt.axhline(cumulative[-1], ls='-.', color='green', label="Total sum")
                plt.axhline(cumulative[-1]*threshold, ls='-.', color='red', label=f'Total sum $\\times$ {threshold}')
                plt.axvline(number_cutoff, ls='--', color='red', 
                            label=f'Cutoff at\nnumber = {number_cutoff:.2f};\n{threshold_new*100:.2f}% reached')
                plt.legend(loc=4)
                plt.subplot(121)
                plt.axvline(number_cutoff, ls='--', color='red')
                #plt.tight_layout()
                plt.subplots_adjust()                 
                plt.show()
                        

            # making masks and applying them to images
            
            nmhg_radial = nmhg
            nmhg1_reserve = nmhg1
            
            filter_mask = nmhg <= number_cutoff
            nmhg = nmhg*filter_mask
            
            filter_mask1 = nmhg1 <= number_cutoff
            nmhg1 = nmhg1*filter_mask1
            
            #plt.imshow(np.rot90(filter_mask1))
            #plt.show()            
            
            if draw_additional:
            
                plt.figure(figsize=(11,5))
                plt.subplot(121)
                plt.title("nmhg1 (filtered)")
                plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.subplot(122)
                plt.title("nmhg (filtered)")
                plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.show()
            
            dddfff["RA_pix"] =  (dddfff["RA"]  - cntr[0] + R)*3600/ang_res # - 1
            dddfff["DEC_pix"] = (dddfff["DEC"] - cntr[1] + R)*3600/ang_res # - 1
            
            dddfff["RA_pix"] = dddfff["RA_pix"].astype(int)
            dddfff["DEC_pix"] = dddfff["DEC_pix"].astype(int)
            
            dddfff["stay"] = filter_mask1[dddfff["RA_pix"], dddfff["DEC_pix"]]
            dddfff = dddfff[dddfff["stay"] == True]
                        
            if draw_additional:
            
                plt.figure(figsize=(11,5))
                plt.subplot(121)
                NMHG2, NMHG_X, NMHG_Y , _ = plt.hist2d(dddfff["RA"], dddfff["DEC"],
                                             bins=len(nmhg1),
                                             norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
                axesforsmooth = list(plt.gca().get_xlim()) + list(plt.gca().get_ylim())
                #print(axesforsmooth)
                plt.gca().set_aspect('equal', 'box')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.title('nmhg1 but in RA/DEC')
                plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
                                                
            dddfff = dddfff.drop("stay", axis=1) 
            dddfff = dddfff.drop("RA_pix", axis=1)
            dddfff = dddfff.drop("DEC_pix", axis=1)
            
            # number_of_unfiltered_photons
            
            number_of_filtered_photons = len(dddfff)
            percent_of_photons = number_of_filtered_photons/number_of_unfiltered_photons
                                
            if draw_additional and True:     
                plt.show()
                
                NMHG2 = np.rot90(NMHG2)
                plt.imshow(convolve(NMHG2, fltr), 
                           norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                           origin='upper',
                           extent = axesforsmooth)
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.gca().set_aspect('equal', 'box')
                plt.title('nmhg1 but in RA/DEC convolved')
                plt.show()
                   
    # this goes to final panels image centered at cntr (which is set to centroid position as default)
           
    SLICE1["todraw"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size),
                                    True, False)
    SLICE1 = SLICE1[SLICE1['todraw'] == True]
    SLICE1 = SLICE1.drop("todraw", axis=1)
        
    SLICE1 = SLICE1[SLICE1['ENERGY']>0.3]      
    SLICE1 = SLICE1[SLICE1['ENERGY']<2.3]
        
    if delete_bright_regions:
                
        SLICE1["RA_pix"] = ((SLICE1["RA"] - cntr[0] + half_size)*3600/ang_res).astype(int) - 1
        SLICE1["DEC_pix"] = ((SLICE1["DEC"] - cntr[1] + half_size)*3600/ang_res).astype(int) - 1
    
        SLICE1['todraw'] = filter_mask[SLICE1["RA_pix"], SLICE1["DEC_pix"]]
        SLICE1 = SLICE1[SLICE1['todraw'] == True]
        
    nmhg, nmhg_x, nmhg_y = np.histogram2d(SLICE1["RA"], SLICE1["DEC"],
                                          bins=int(2*half_size*3600/ang_res),
                                          #norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                          range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                          (cntr[1]-half_size, cntr[1]+half_size)]))
                                                            
    axesforsmooth = [nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]] 
   
    if draw:
           
        trtr = plt.imshow(convolve(np.rot90(nmhg), Gaussian2DKernel(1)), 
                          norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                          origin='upper',
                          extent = axesforsmooth)     
        
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue')
        plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
        
        if delete_superfluous:
                                                        
            for vv in vicenter:
                plt.scatter(vv[0], vv[1], color='red', label = 'Subhaloes', s=7)
                        
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        
        x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
        y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
        y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]   
        plt.plot((x_s+5/60, x_s-5/60), (y_s, y_s), color='white')
        plt.text(x_s, y_S, f'10 arcmin $\\approx$ {10/60*D_A*np.pi/180:.0f} kpc', 
                 color='white', ha='center', va='center')
        
        plt.xlim(cntr[0]-half_size, cntr[0]+half_size)
        plt.ylim(cntr[1]-half_size, cntr[1]+half_size)
        plt.gca().set_aspect('equal', 'box')
        
        plt.xlabel("RA, deg", size=13)
        plt.ylabel("DEC, deg", size=13)
        plt.xticks(size=13)
        plt.yticks(size=13)
        cb = plt.colorbar(trtr, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=13)
        cb.set_label(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin", size=13)
        ttiittllee = f'#{current_cluster_number}: z={ztrue:.3f}, A={AREA:.1f} min$^2$'
        
        if not delete_bright_regions:
            plt.title(ttiittllee, fontsize=15)
            #cb.ax.set_yticklabels(['0', '1', '10', '100'])
        else:
            plt.title(ttiittllee+f', RP={100*percent_of_photons:.1f}%', fontsize=15)
            #cb.ax.set_yticklabels(['0', '1'])
        
        handles, labels = plt.gca().get_legend_handles_labels()
        #l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label=str(r)+"$\cdot R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l2])
        #plt.legend(handles=handles, loc=3, fontsize=13)
        #plt.show()
              
    #if redshifted_back:
    #    return dddfff.mul(1+ztrue)
    #else:
    #    return dddfff

    return nmhg
 
def kruzhok(r_pixels, mm, NMHG, d_pixels):

    kusok = np.zeros((2*d_pixels+1, 2*d_pixels+1))
        
    for i in range(mm[0]-d_pixels, mm[0]+d_pixels+1):
        for j in range(mm[1]-d_pixels, mm[1]+d_pixels+1):
            kusok[i - (mm[0]-d_pixels)][j - (mm[1]-d_pixels)] = NMHG[i][j]
        
    Y, X = np.ogrid[(mm[1]-d_pixels):(mm[1]+d_pixels+1), (mm[0]-d_pixels):(mm[0]+d_pixels+1)]
    dist_from_center = np.sqrt((X - mm[0])**2 + (Y-mm[1])**2)
    
    mask = dist_from_center <= r_pixels
        
    return mask*kusok, mask
    

def draw_84_panels(del_br_reg):

    NNN = 84
    
    #size = 6

    #plt.figure(figsize=((size)*7+6*3, 5*12+11*2.5))
    plt.tight_layout()
    
    for cl_num in tqdm(clusters.index[:NNN]):
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        pho_list = extract_photons_from_cluster(cl_num, r = 1, draw=True, delete_bright_regions=del_br_reg)


def calc_l_T(T, T_left, T_right, Xplot=False):
  
    x.Xset.chatter = 0
    x.AllData.clear()
    x.AllData.removeDummyrsp()
    x.AllData.dummyrsp(lowE=0.1, highE=10.0, nBins=1024)
    x.Xset.addModelString("APEC_TRACE_ABUND", "0")

    if Xplot:
        x.Plot.device = '/xs'
    else:
        x.Plot.device = '/null'
    
    expo = 40000
    Ab = 0.3
    Norm = 1
    z = 0.0
    nH = 0.01
    
    mod = x.Model('phabs*apec')
    mod.setPars(nH, T, Ab, z, Norm)
    
    x.Plot("model")
    x.Plot.xAxis = "keV"
    x.AllModels.show()
    
    x.AllModels.calcFlux(f"{T_left} {T_right}")
    flx1 = x.AllModels(1).flux[0]   # ergs/cm^2
    flx2 = x.AllModels(1).flux[3]   # photons 
    
    RMF_NAME = '../erosita/erosita_pirmf_v20210719.rmf'
    ARF_NAME = '../erosita/tm1_arf_open_000101v02.fits'
    
    fs = x.FakeitSettings(response = RMF_NAME, 
                               arf = ARF_NAME, 
                        background = '', 
                          exposure = expo, 
                        correction = '', 
                      backExposure = '', 
                          fileName = 'fakeit.pha')
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = True)

    x.AllData.ignore(f"**-{T_left} {T_right}-**")             # IMPORTANT !
    x.AllData.show()
    #x.AllModels.setEnergies("reset")
    
    x.Plot("data")
    x.Plot.xAxis = "keV"
    #xVals = x.Plot.x()
    #yVals = x.Plot.y()
    
    cr = x.AllData(1).rate[2]
    
#    ens = x.AllData(1).energies
#    E_i = np.zeros(len(ens))
#    dE = np.zeros(len(ens))   
#    for i in ens:
#        dE[ens.index(i)] = i[1]-i[0]
#        E_i[ens.index(i)] = (i[0]+i[1])/2    
#    s_i = x.AllData(1).values
    
    return flx1, flx2, cr #, np.dot(xVals, yVals)/sum(yVals), np.dot(E_i, s_i)/cr
