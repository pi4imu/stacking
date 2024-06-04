# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r, centroid=True, delete_bright_regions=True, draw=True, draw_additional=False, redshifted_back=True):

    # there are several cases of SAME ihal for DIFFERENT cluster numbers
    # this is the reason for using cluster number as a counter
    
    current_cluster = clusters.loc[current_cluster_number]
    
    RA_c = current_cluster["x_pix"]*30-5
    DEC_c = current_cluster["y_pix"]*30-5
    R_vir = current_cluster["Rrel"]*30
    R_500 = current_cluster["R500"]*0.704  # kpc
    ztrue = current_cluster["z_true"]
    
    D_A = 343000 # FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    R_500_rescaled = R_500/D_A*180/np.pi *2 # degrees
    
    snap_id_str = binned_clusters[current_cluster_number][1]   # id of photon list
        
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
    
    SLICE = t.to_pandas()        # for photons extraction
    SLICE1 = t.to_pandas()       # for drawing
    SLICE2 = t.to_pandas()       # for center searching if it is not from table
    SLICE3 = t.to_pandas()       # for rescaling SLICE2
    
    if r == 'Rvir':
        R = R_vir
        print('Really?')
    elif r == 'R500':
        R = R_500_rescaled
        print('OBSOLETE! Change R500 to numerical value!')
    else:
        R = r * R_500_rescaled
        
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
        half_size = 5*R_500_rescaled
        
        if (current_cluster_number != 13334) and (current_cluster_number != 18589):
            hs4s = half_size/3/5
        else:
            hs4s = half_size/1
        # making 2D histogram of size 2*half_size with center (RA_c, DEC_c) without drawing
                
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < hs4s) & (np.abs(SLICE2["DEC"]-DEC_c) < hs4s), True, False)
        whattodraw = SLICE2[SLICE2['what'] == True]
        whattodraw = whattodraw.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], bins=int(2*hs4s*3600/ang_res))
                       
        # centroid position
        
        #nmhg = np.zeros((100,100))
        #for i in range(0, len(nmhg)):
        #    for j in range(0, len(nmhg)):
        #        if i>=70:
        #            nmhg[i,j] = 1
        #print(nmhg)
        
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
        
        # searching for the coordinates of maximum value
        #whereee = np.concatenate(np.where(nmhg == max(nmhg.flatten())))
        #reeeversed = [a*ang_res/60/60 for a in whereee]
        #xeeec = plt.gca().get_xlim()[0] + reeeversed[0]
        #yeeec = plt.gca().get_ylim()[0] + reeeversed[1]
        #print(xeeec)
        #print(yeeec)
        #print(np.where(nmhg == max(nmhg.flatten())))
        #print(max(nmhg.flatten()))
        #m = [int(a[0]) for a in np.where(nmhg == max(nmhg.flatten()))]
        #plt.scatter(xeeec, yeeec, color='dodgerblue', label = 'Max value')
        #plt.gca().add_patch(plt.Circle((xeeec, yeeec), R_500_rescaled, color='yellow', 
        #linestyle="--", lw=3, fill = False))                
        
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
            
            #c2 = np.array( ) + shift #+ shift2
            
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
            
            threshold = 0.95
            
            beeens = np.geomspace(1, np.max(nmhg1.flatten()), 50)
            amount_in_bin, bin_borders = np.histogram(nmhg1.flatten(), bins = beeens)
            #bin_centers = [int(b) for b in bin_borders[:-1]]
            bin_centers = (bin_borders[:-1]+bin_borders[1:])/2
            #print(bin_centers)
            cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))]
            number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
            index_cutoff = list(bin_centers).index(number_cutoff)
            threshold_new = sum(amount_in_bin[:index_cutoff]) / cumulative[-1]
            
            #sum_amount, i = 0, 0
            #while sum_amount <= threshold * sum(amount_in_bin*bin_centers):
            #    sum_amount = sum_amount + (amount_in_bin*bin_centers)[i]
            #    #print(i, bin_centers[i], sum_amount/sum(amount_in_bin*bin_centers))
            #    i = i + 1
            #threshold = (sum_amount-(amount_in_bin*bin_centers)[i-1])/sum(amount_in_bin*bin_centers)        
            #number_cutoff = bin_centers[i-1]

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
            
            # another method of excluding bright pixels (by making rings):
            
            if False:
                
                ffltr =  nmhg1_reserve*0 #(nmhg1_reserve == 0)
                
                #plt.imshow(ffltr)
                #plt.show()
                
                r_pixels_max = int(R*3600/ang_res)
                brightness = []   #[k0[0].sum()/sum(k0[1].flatten())]
            
                ring_width = 20
                                  
                for rr in range(0, r_pixels_max+1, ring_width):
        
                    k1 = kruzhok(rr, c1, nmhg_radial, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c1, nmhg_radial, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    koltso = (sum(k2[1].flatten())-sum(k1[1].flatten()))
                    print('3.14 * (', (rr + ring_width)**2, '-', rr**2, ') =', koltso)
                     
                    plt.figure(figsize=(19,3))                
                    plt.subplot(151)
                    ring = np.rot90(ring)
                    plt.imshow(ring, 
                               norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                               origin='upper')
                    plt.colorbar(fraction=0.046, pad=0.04)
                    
                    plt.subplot(152)
                    bbiinnss = np.geomspace(1, np.max(ring.flatten()), 50)
                    bbiinnss = sorted(bbiinnss)
                    #print(ring.flatten(), bbiinnss)
                    amount_in_bin, bin_borders = np.histogram(ring.flatten(), bins = bbiinnss)
                    #print(amount_in_bin, bin_borders)                    
                    amount_in_bin = amount_in_bin/koltso
                    bin_centers = (bin_borders[:-1]+bin_borders[1:])/2
                    #bin_centers = [int(b) for b in bin_borders[:-1]]
                    plt.plot(bin_centers, amount_in_bin)
                    
                    plt.subplot(153)
                    
                    cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))]
                    number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
                    index_cutoff = list(bin_centers).index(number_cutoff)
                    threshold_new = sum(amount_in_bin[:index_cutoff]) / cumulative[-1]
                    
                    plt.plot(bin_centers, cumulative)
                    plt.axhline(cumulative[-1], ls='-.', color='green', label="Total sum")
                    plt.axhline(cumulative[-1]*threshold, ls='-.', color='red', label=f'Total sum $\\times$ {threshold}')
                    plt.axvline(number_cutoff, ls='--', color='red', 
                                label=f'Cutoff at\nnumber = {number_cutoff:.2f};\n{threshold_new*100:.2f}% reached')
                    plt.legend(loc=4)
                    
                    plt.subplot(152)
                    plt.axvline(number_cutoff, ls='--', color='red')
                    
                    plt.subplot(154)
                    filter_mask2 = ring <= number_cutoff
                    ring = ring*filter_mask2
                    filter_mask3 = filter_mask2 == 0
                    plt.imshow(filter_mask3)
                    
                    plt.subplot(155)
                    ffltr = ffltr + ring

                    plt.imshow(ffltr, 
                               norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                               origin='upper')
                    plt.colorbar(fraction=0.046, pad=0.04)
                    #plt.subplots_adjust()
                    plt.tight_layout()
                    plt.show()
            
                nmhg1 = np.rot90(ffltr, 3)
            
            
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
               
            # here we obtain brightness profile inside R_500_rescaled               
               
            if False:
            
                plt.subplot(122)
              
                r_pixels_max = int(R*3600/ang_res)
               
                #k0 = kruzhok(0, c, nmhg, r_pixels_max+1)
                brightness = []   #[k0[0].sum()/sum(k0[1].flatten())]
                brightness_filtered = []
            
                ring_width = 10
                 
                for rr in range(0, r_pixels_max+1):
        
                    k1 = kruzhok(rr, c1, nmhg, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c1, nmhg, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    brightness.append(ring.sum()/(sum(k2[1].flatten())-sum(k1[1].flatten())))
                    
                    k1 = kruzhok(rr, c, nmhg_radial, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c, nmhg_radial, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    brightness_filtered.append(ring.sum()/(sum(k2[1].flatten())-sum(k1[1].flatten())))                    
            
                r500_pix = int(R_500_rescaled*3600/ang_res)                
                plt.plot(np.linspace(0, r_pixels_max+1, r_pixels_max+1)/r500_pix, brightness)
                plt.plot(np.linspace(0, r_pixels_max+1, r_pixels_max+1)/r500_pix, brightness_filtered)
                #plt.axvline((brightness.index(max(brightness))+1)/r500_pix, ls='--', color='black')
                plt.xlabel("Radius in units of $R_{500}$")
                plt.ylabel("Brightness in relative units")
                #plt.axhline(brightness_max, ls='--', color='red', 
                #label=f'{threshold*100:.2f} % cutoff\nat brightness = {brightness_max:.2f}')
                plt.xscale("log")
                plt.yscale("log")
                #plt.legend()
                plt.subplots_adjust()                 
                plt.tight_layout()
                                
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
    
    if draw:
               
        SLICE1["todraw"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size),
                                    True, False)
        SLICE1 = SLICE1[SLICE1['todraw'] == True]
        SLICE1 = SLICE1.drop("todraw", axis=1)
        
        #SLICE1 = SLICE1[SLICE1['ENERGY']>1.3]      
        #№SLICE1 = SLICE1[SLICE1['ENERGY']<2.3]
        
        if delete_bright_regions:
                
            SLICE1["RA_pix"] = ((SLICE1["RA"] - cntr[0] + half_size)*3600/ang_res).astype(int) - 1
            SLICE1["DEC_pix"] = ((SLICE1["DEC"] - cntr[1] + half_size)*3600/ang_res).astype(int) - 1
        
            #plt.hist2d(SLICE1["RA_pix"], SLICE1["DEC_pix"], 
            #           bins=len(nmhg),
            #           norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
            #plt.gca().set_aspect('equal', 'box')
            #plt.colorbar(fraction=0.046, pad=0.04)
            #plt.show()
        
            #print(min(SLICE1["RA_pix"]), min(SLICE1["DEC_pix"]))
        
            SLICE1['todraw'] = filter_mask[SLICE1["RA_pix"], SLICE1["DEC_pix"]]
            SLICE1 = SLICE1[SLICE1['todraw'] == True]
        
        if False:
        
            nmhg, _, _, trtr = plt.hist2d(SLICE1["RA"], SLICE1["DEC"],
                                          bins=int(2*half_size*3600/ang_res),
                                          norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                          range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                          (cntr[1]-half_size, cntr[1]+half_size)]))
        
        else:
        
            nmhg, nmhg_x, nmhg_y = np.histogram2d(SLICE1["RA"], SLICE1["DEC"],
                                              bins=int(2*half_size*3600/ang_res),
                                              #norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                              range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                              (cntr[1]-half_size, cntr[1]+half_size)]))
            
            #print([nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]])                                                  
            axesforsmooth = [nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]] 
           
            trtr = plt.imshow(convolve(np.rot90(nmhg), Gaussian2DKernel(1)), 
                              norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                              origin='upper',
                              extent = axesforsmooth)
                        
        
        # obsolete (it was needed for estimation of correct position of circles)
        
        #m_x, m_y = RA_c - half_size + c[0]*ang_res/3600, DEC_c - half_size + c[1]*ang_res/3600
        #r_degrees = r_pixels_max*ang_res/3600
        #plt.gca().add_patch(plt.Circle((m_x, m_y), r_degrees, color='red', 
        #linestyle="-", lw=1, fill = False, label = 'Brightest'))
        #plt.gca().add_patch(plt.Rectangle((m_x-r_degrees, m_y-r_degrees), 
        #2*r_degrees, 2*r_degrees, color='red', linestyle="-", lw=1, fill = False, label = 'Brightest'))
                    
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue')
        plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
            
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        
        x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
        y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
        y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]
        #plt.scatter(x_s, y_s, color='red')     
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
            cb.ax.set_yticklabels(['0', '1', '10', '100'])
        else:
            #plt.title(ttiittllee+f'\nPercentage of remaining photons: {100*percent_of_photons:.1f}%', fontsize=15)
            plt.title(ttiittllee+f', RP={100*percent_of_photons:.1f}%', fontsize=15)
            cb.ax.set_yticklabels(['0', '1'])
        
        handles, labels = plt.gca().get_legend_handles_labels()
        #l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label=str(r)+"$\cdot R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l2])
        plt.legend(handles=handles, loc=3, fontsize=13)
        #plt.show()
                
              
    if redshifted_back:
        return dddfff.mul(1+ztrue)
    else:
        return dddfff

 
def kruzhok(r_pixels, mm, NMHG, d_pixels):

    kusok = np.zeros((2*d_pixels+1, 2*d_pixels+1))
        
    for i in range(mm[0]-d_pixels, mm[0]+d_pixels+1):
        for j in range(mm[1]-d_pixels, mm[1]+d_pixels+1):
            kusok[i - (mm[0]-d_pixels)][j - (mm[1]-d_pixels)] = NMHG[i][j]
        
    Y, X = np.ogrid[(mm[1]-d_pixels):(mm[1]+d_pixels+1), (mm[0]-d_pixels):(mm[0]+d_pixels+1)]
    dist_from_center = np.sqrt((X - mm[0])**2 + (Y-mm[1])**2)
    
    mask = dist_from_center <= r_pixels
        
    return mask*kusok, mask
    

def draw_84_panels(mode):

    NNN = 84
    
    if mode=='IMAGE':
        size = 6
    else:
        size = 5

    plt.figure(figsize=((size)*7+6*3, 5*12+11*2.5))
    #plt.figure(figsize=((size)*3+6*3, (size)*4+11*2.5))
    plt.tight_layout()
    
    #if mode!='IMAGE':
    
        #temp_compare = {}
        #lumin_compare = {}
        #average_ene = {}
        
    for cl_num in clusters.index[:NNN]:
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        if mode=='IMAGE':
        
            pho_list = extract_photons_from_cluster(cl_num, r = 1, draw=True, delete_bright_regions=True)
        
        else:
        
            #cl_T500 = clusters.loc[cl_num]["T500"]
            #cl_lum = clusters.loc[cl_num]["Lx500"]
    
            SP = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=False, inside_radius="R500",
                                            dbr=False, Xplot=False, plot=True, draw_only=mode)

            #temp_compare[cl_num] = [cl_T500, SP[0][:3]]
            #lumin_compare[cl_num] = [cl_lum, SP[1][:3]]
            #average_ene[cl_num] = [SP[2]]
