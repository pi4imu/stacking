# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r=1.0, centroid=True, delete_superfluous=False, draw=True, histlen=2001, withagn=False, ARF_weights=True):

    # there are several cases of SAME ihal for DIFFERENT cluster numbers
    # this is the reason for using cluster number as a counter
    
    current_cluster = clusters.loc[current_cluster_number]
    
    snap_id_str = binned_clusters[current_cluster_number][1]   # id of photon list

    if snap_id_str == '124':
        zslice = 0.174192889973847
    elif snap_id_str == '128':
        zslice = 0.13708140389145
    elif snap_id_str == '132':
        zslice = 0.101142861718869
    elif snap_id_str == '136':
        zslice = 0.0663401914452304
    elif snap_id_str == '140':
        zslice = 0.032637492755919
    
    RA_c = current_cluster["x_pix"]*30-5
    DEC_c = current_cluster["y_pix"]*30-5
    R_vir = current_cluster["Rrel"]*30
    R_500 = current_cluster["R500"] / 0.704 / (1+zslice) # kpc
    ztrue = current_cluster["z_true"]
    
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    R_500_rescaled = R_500/D_A.value*180/np.pi # degrees
    
    # angularDiameterDistance( z, H0, omegaM, omegaLambda ) angularDiameterDistance( z, H0, omegaM, omegaLambda )dad angularDiameterDistance( z, H0, omegaM, omegaLambda ) 
    
    # parseFloat(col10) / (1+parseFloat(col7)) / 0.704 / angularDiameterDistance( parseFloat(col7), 70.4, 0.272, 0.728 ) / 1000 * 180 / 3.141592
            
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
   
   # adding the background from AGN only (obsolete and needs to be optimized)
    
   # t1 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_124.fits", hdu=2)
   # t2 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_128.fits", hdu=2)
   # t3 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_132.fits", hdu=2)
   # t4 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_136.fits", hdu=2)
   # t5 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_140.fits", hdu=2)    
   # t = vstack([t1, t2, t3, t4, t5])
   # if withagn:
   #     for snapid in [124, 128, 132, 136, 140]:
   #         t1 = Table.read("../data/eROSITA_30.0x30.0/Phox/AGNphlist_"+str(snapid)+".fits", hdu=2)
   #         t = vstack([t, t1])
    
    SLICE = t.to_pandas()        # for photons extraction
    SLICE1 = t.to_pandas()       # for drawing
    SLICE2 = t.to_pandas()       # for center searching if it is not from table
    
    R = r * R_500_rescaled    
    AREA = np.pi*R**2*3600   # min2
    
    #setting area and resolution for searching for center
    
    ang_res = 4                    # this value will be changed later 
    halfsidelength = 10             # in units of R500
    half_size = halfsidelength*R   # in degrees
            
    if not centroid:
    
        # taking photons from circle centered at RA_c, DEC_c
    
        SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        cntr = (RA_c, DEC_c) # for drawing
    
    else:
    
        # hs4s - half size for search
            
        if (current_cluster_number != 13334) and (current_cluster_number != 18589):
            hs4s = R
        else:
            hs4s = 3*R
        # for peak pixel    
        if (current_cluster_number in [820, 1819, 6740, 10551, 11272, 13675, 17017, 18589, 18610]):
            hs4s = R/4    
        
        # making 2D histogram with side length 2*half_size with center (RA_c, DEC_c) without drawing
                
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < hs4s) & (np.abs(SLICE2["DEC"]-DEC_c) < hs4s), True, False)
        search4centroid = SLICE2[SLICE2['what'] == True]
        search4centroid = search4centroid.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(search4centroid["RA"], search4centroid["DEC"], bins=int(2*hs4s*3600/ang_res))
        
        smooth = 40 / 3600 / R_500_rescaled * 100
        #print("R_500 =", R_500_rescaled, "degrees;   kernel =", smooth, "pixels")
        nmhg = convolve_fft(nmhg, Gaussian2DKernel(smooth))
                       
        # centroid position
                
        psum = sum(nmhg.flatten())
        c_x, c_y = 0, 0
        
        for i in range(0, len(nmhg)):
            for j in range(0, len(nmhg)):
                c_x = c_x + i*nmhg[i,j]/psum
                c_y = c_y + j*nmhg[i,j]/psum
        
        # position of centroid in units of pixels relative to the upper left border    
        c = [int(c_x), len(nmhg)-int(c_y)]
                        
        c_x_1 = RA_c - hs4s + c_x*ang_res/3600 
        c_y_1 = DEC_c - hs4s + c_y*ang_res/3600
        cntr = (c_x_1, c_y_1) # position in degrees
        
        # taking photons from circle centered at centroid
        
        #SLICE["check"]=np.where((SLICE["RA"]-c_x_1)**2+(SLICE["DEC"]-c_y_1)**2 <= R**2, True, False)
        #df = SLICE[SLICE['check'] == True]
        #dddfff = df.drop("check", axis=1)
    
        # searching for the coordinates of maximum value
        
        if True:
            #print(max(nmhg.flatten()))
            maxpixpos = np.where(nmhg == max(nmhg.flatten()))
            #print(maxpixpos)
            #print(nmhg[maxpixpos])
            whereee = [maxpixpos[0][0], maxpixpos[1][0]]  # np.concatenate(maxpixpos)
            #print(whereee)
            #print(nmhg[whereee[0], whereee[1]])
            xeeec = RA_c  - hs4s + whereee[0]*ang_res/3600
            yeeec = DEC_c - hs4s + whereee[1]*ang_res/3600
            #print(xeeec, yeeec)
            cntr = (xeeec, yeeec)

            #plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=0.1, linscale=1), origin='upper')
            #plt.scatter(whereee[0], len(nmhg)-whereee[1], c='r')
            #plt.scatter(whereee[3], whereee[1], c='r')
            #plt.colorbar(fraction=0.046, pad=0.04)
            #plt.show()

            SLICE["check"]=np.where((SLICE["RA"]-xeeec)**2+(SLICE["DEC"]-yeeec)**2 <= R**2, True, False)
            df = SLICE[SLICE['check'] == True]
            dddfff = df.drop("check", axis=1)                   
    
    # this goes to final panels image centered at cntr (which is set to centroid position by default)
           
    SLICE1["todraw"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size),
                                    True, False)
    SLICE1 = SLICE1[SLICE1['todraw'] == True]
    SLICE1 = SLICE1.drop("todraw", axis=1)
    
    # setting the energy band
        
    SLICE1 = SLICE1[SLICE1['ENERGY']>0.3]           # previosly 0.3-2.3
    SLICE1 = SLICE1[SLICE1['ENERGY']<2.3]
    
    # making histogram without or with applying weights that are given by ARF
    
    if not ARF_weights:
     
        nmhg, nmhg_x, nmhg_y = np.histogram2d(SLICE1["RA"], SLICE1["DEC"],
                                              bins=histlen, #int(2*half_size*3600/ang_res),
                                              #norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                              range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                              (cntr[1]-half_size, cntr[1]+half_size)])) #, density=True)
                                              #weights=SLICE1["ENERGY"]
        LINTHRESH = 0.000001                                      
        
    else:
    
        arf_survey = fits.open('../erosita/esf10.Dsur1234regR3cCaXv2.0001.arf')[1].data
        
        sl = SLICE1
        #sl["FLUX"] = sl["ENERGY"] / 1000 / 10000 / 4**2 * 60**2      # keV/cm2/s/arc
        sl["EFF_AREA"] = np.interp(sl["ENERGY"], arf_survey["ENERG_LO"], arf_survey["SPECRESP"]) # cm2
        #sl["RATE"] = sl["FLUX"] * sl["EFF_AREA"]    # keV/s
        
        nmhg, nmhg_x, nmhg_y = np.histogram2d(sl["RA"], sl["DEC"],
                                              bins=histlen, #int(2*half_size*3600/ang_res),
                                              #norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                              range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                              (cntr[1]-half_size, cntr[1]+half_size)]),
                                              weights=sl["EFF_AREA"], density=False)
        LINTHRESH = 0.0001
        
    axesforsmooth = [nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]]    # for drawing
      
    ang_res = 2*half_size*3600 / histlen               # new value instead of 4'' by default    

    # important rescaling
    
    nmhg = nmhg / 1000 / 10000 / ang_res**2 * 60**2
    
    f1 = 1000/R_500
    f2 = E(ztrue)**(-4)*(1+ztrue)**3
    factor = f1*f2
    nmhg = nmhg*factor
        
    # deleting haloes in the vicinity
        
    if delete_superfluous:
        
        # taking haloes only from the slice to which this cluster belongs 
        
        VICINITY = np.where( ((clusters_all["x_pix"]*30-5 - c_x_1)**2 + 
                              (clusters_all["y_pix"]*30-5 - c_y_1)**2 < 2*half_size**2) ) # & 
                             # ( np.abs(clusters_all["z_true"] - ztrue) < 0.027141 ) )  # 1 changed to 2
        #print(VICINITY)
        vclu = clusters_all.loc[VICINITY]  
        #display(vclu)
                
        for clcl in VICINITY:
        
            vicinity_current = clusters_all.loc[clcl].drop(current_cluster_number)

            D_A_ = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(vicinity_current["z_true"])*1000 # kpc
            R_500_ = vicinity_current["R500"]/ 0.704 / (1+zslice)
            R_500_s = R_500_/D_A_*180/np.pi # degrees
            
            vicenter = list(zip(vicinity_current["x_pix"].values*30-5, 
                                vicinity_current["y_pix"].values*30-5, 
                                R_500_s))
            #print(clcl)
            #print(vicenter)
   
        # creating the mask which will go to output along with unfiltered field
        
        nmhg_mask = np.zeros((histlen, histlen))
        #MASK = pd.DataFrame([])
        
        for vic in vicenter:
        
            #print(vic)
            ra_pix = np.abs(vic[0]-cntr[0]-half_size)*3600/ang_res
            dec_pix = (vic[1]-cntr[1]+half_size)*3600/ang_res
            radius_pix = vic[2]*3600/ang_res
            
            ra_pix = ra_pix.astype(int)
            dec_pix = dec_pix.astype(int)
            radius_pix = radius_pix.astype(int)
            
            circle_mask = create_circle_mask(ra_pix, dec_pix, 1.6*radius_pix, len(nmhg_mask))  # R200c
            #if circle_mask[int(histlen/2), int(histlen/2)] == 0:
            nmhg_mask = nmhg_mask + circle_mask
        
        # manual filtering

        if True:
        
            # something big and not very near
            
            if (current_cluster_number == 17638):
                pup = create_circle_mask(550, 950, 200, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-550)*ang_res/3600-half_size+cntr[0], 
                                 950*ang_res/3600-half_size+cntr[1], 
                                 200/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup
                pup = create_circle_mask(850, 1200, 150, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-850)*ang_res/3600-half_size+cntr[0], 
                                 1200*ang_res/3600-half_size+cntr[1], 
                                 150/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup
                            
            if (current_cluster_number == 18589):
                pup = create_circle_mask(700, 1080, 80, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-700)*ang_res/3600-half_size+cntr[0], 
                                 1080*ang_res/3600-half_size+cntr[1], 
                                 80/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup
                
            if (current_cluster_number == 10018):
                pup = create_circle_mask(1350, 1070, 100, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-1350)*ang_res/3600-half_size+cntr[0], 
                                 1070*ang_res/3600-half_size+cntr[1], 
                                 100/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup
                
            if (current_cluster_number == 13675):
                pup = create_circle_mask(1500, 1050, 100, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-1500)*ang_res/3600-half_size+cntr[0], 
                                 1050*ang_res/3600-half_size+cntr[1], 
                                 100/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup        

# because filtering by haloes in all slices has deleted these clumps already
                
#            if (current_cluster_number == 17421):
#                pup = create_circle_mask(450, 1700, 70, 2001)
#                vicenter.append(((2000-450)*ang_res/3600-half_size+cntr[0], 
#                                 1700*ang_res/3600-half_size+cntr[1], 
#                                 70*ang_res/3600))
#                nmhg_mask = nmhg_mask + pup
                
            if False:  # nearest
                            
                if (current_cluster_number == 7996):
                    pup = create_circle_mask(1050, 1200, 100, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-1050)*ang_res/3600-half_size+cntr[0], 
                                     1200*ang_res/3600-half_size+cntr[1], 
                                     100/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup

                if (current_cluster_number == 14857):
                    pup = create_circle_mask(880, 1150, 100, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-880)*ang_res/3600-half_size+cntr[0], 
                                     1150*ang_res/3600-half_size+cntr[1], 
                                     100/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup

                if (current_cluster_number == 1819):
                    pup = create_circle_mask(1200, 1050, 100, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-1200)*ang_res/3600-half_size+cntr[0], 
                                     1050*ang_res/3600-half_size+cntr[1], 
                                     100/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup                   
    
                if (current_cluster_number == 11141):
                    pup = create_circle_mask(850, 1200, 120, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-850)*ang_res/3600-half_size+cntr[0], 
                                     1200*ang_res/3600-half_size+cntr[1], 
                                     120/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup  

                if (current_cluster_number == 7308):
                    pup = create_circle_mask(870, 870, 70, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-870)*ang_res/3600-half_size+cntr[0], 
                                     870*ang_res/3600-half_size+cntr[1], 
                                     70/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup 

                if (current_cluster_number == 9836):
                    pup = create_circle_mask(1100, 1150, 100, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-1100)*ang_res/3600-half_size+cntr[0], 
                                      1150*ang_res/3600-half_size+cntr[1], 
                                      100/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup                 

                if (current_cluster_number == 9240):
                    pup = create_circle_mask(1050, 1150, 100, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-1050)*ang_res/3600-half_size+cntr[0], 
                                      1150*ang_res/3600-half_size+cntr[1], 
                                      100/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup 

                if (current_cluster_number == 17086):
                    pup = create_circle_mask(900, 1100, 70, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-900)*ang_res/3600-half_size+cntr[0], 
                                      1100*ang_res/3600-half_size+cntr[1], 
                                      70/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup
    
                    pup = create_circle_mask(800, 1250, 70, 2001)
                    pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                    vicenter.append(((2000-800)*ang_res/3600-half_size+cntr[0], 
                                      1250*ang_res/3600-half_size+cntr[1], 
                                      70/1.6*ang_res/3600))
                    nmhg_mask = nmhg_mask + pup                 
            
#            if (current_cluster_number == 4613):
#                pup = create_circle_mask(380, 1350, 70, 2001)
#                vicenter.append(((2000-380)*ang_res/3600-half_size+cntr[0], 
#                                 1350*ang_res/3600-half_size+cntr[1], 
#                                 70*ang_res/3600))
#                nmhg_mask = nmhg_mask + pup              
            
            # periphery
            
            if (current_cluster_number == 17421):
                pup = create_circle_mask(80, 820, 70, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-80)*ang_res/3600-half_size+cntr[0], 
                                 820*ang_res/3600-half_size+cntr[1], 
                                 70*ang_res/3600))
                nmhg_mask = nmhg_mask + pup  

            if (current_cluster_number == 4613):
                pup = create_circle_mask(920, 1940, 70, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-920)*ang_res/3600-half_size+cntr[0], 
                                 1940*ang_res/3600-half_size+cntr[1], 
                                 70/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup  

           # if (current_cluster_number == 171):
           #      pup = create_circle_mask(1650, 680, 70, 2001)
           #      vicenter.append(((2000-1650)*ang_res/3600-half_size+cntr[0], 
           #                       680*ang_res/3600-half_size+cntr[1], 
           #                       70*ang_res/3600))
           #      nmhg_mask = nmhg_mask + pup
           
            if (current_cluster_number == 171):
                 pup = create_circle_mask(1600, 1650, 70, 2001)
                 vicenter.append(((2000-1600)*ang_res/3600-half_size+cntr[0], 
                                  1650*ang_res/3600-half_size+cntr[1], 
                                  70/1.6*ang_res/3600))
                 
   #              plt.imshow(pup)
   #              plt.show()
                                                   
                 pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                 
          #       block_size = (old // new, old // new)
                
          #      downsampled_mask = block_reduce(pup_trimmed, block_size=block_size, func=np.any)
          #       pup = downsampled_mask[:new, :new]
                 
                 nmhg_mask = nmhg_mask + pup
                 
    #             plt.imshow(pup)
    #             plt.show()

            if (current_cluster_number == 3155):
                pup = create_circle_mask(1950, 1200, 70, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-1950)*ang_res/3600-half_size+cntr[0], 
                                  1200*ang_res/3600-half_size+cntr[1], 
                                  70/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup     

            if (current_cluster_number == 7191):
                pup = create_circle_mask(930, 1850, 100, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-930)*ang_res/3600-half_size+cntr[0], 
                                  1850*ang_res/3600-half_size+cntr[1], 
                                  100/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup

            if (current_cluster_number == 7553):
                pup = create_circle_mask(1200, 1950, 100, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-1200)*ang_res/3600-half_size+cntr[0], 
                                  1950*ang_res/3600-half_size+cntr[1], 
                                  100/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup

            if (current_cluster_number == 14726):
                pup = create_circle_mask(300, 420, 80, 2001)
                pup = resize(pup.astype(float), (histlen, histlen), 
                              order=3, mode='reflect', anti_aliasing=False, preserve_range=True) > 0.5
                vicenter.append(((2000-300)*ang_res/3600-half_size+cntr[0], 
                                  420*ang_res/3600-half_size+cntr[1], 
                                  80/1.6*ang_res/3600))
                nmhg_mask = nmhg_mask + pup
                                   
     #   kernel = 26 / 3600 / R_500_rescaled * 100
     #   print("R_500 =", R_500_rescaled, "degrees;   kernel =", kernel, "pixels")
        nmhg = convolve_fft(nmhg, Gaussian2DKernel(1))
        # >3 na 201 - uzhe mnogo
        
        nmhg_mask[nmhg_mask > 1] = True   
        nmhg_mask = np.rot90(nmhg_mask)         # some magic
        
        nmhg_unfiltered = nmhg
        nmhg = nmhg_unfiltered * (1-nmhg_mask)
        
        # keeping everything inside R500
        
        center_ones = create_circle_mask(int(histlen/2), int(histlen/2), int(histlen/20), histlen)
        savecenter = nmhg_unfiltered*center_ones
        
        nmhg = nmhg*(1-center_ones) + savecenter
        nmhg_mask = (1-nmhg_mask)*(1-center_ones) + center_ones
        nmhg_mask[nmhg_mask > 1] = True
        nmhg_mask = 1- nmhg_mask
        
    # attempt to take galaxies into account        
            
   #     galaxies_all = pd.read_csv("../data/eROSITA_30.0x30.0/Catalouges/galaxies.dat", sep='\\s+', header=0)    
        VICINITY_GAL = np.where( 
        ((galaxies_all["x_pix"]*30-5 - c_x_1)**2 + (galaxies_all["y_pix"]*30-5 - c_y_1)**2 < 2*half_size**2)
        & ( np.abs(galaxies_all["z_true"] - ztrue) < 0.017141 ) )          
        #vclu_gal = galaxies_all.loc[VICINITY_GAL]           # no need
        
        if False:
            #galaxies_all = pd.read_csv("../data/eROSITA_30.0x30.0/Catalouges/galaxies.dat", sep='\\s+', header=0)    
            VICINITY_GAL = np.where( 
                ((galaxies_all["x_pix"]*30-5 - c_x_1)**2 + (galaxies_all["y_pix"]*30-5 - c_y_1)**2 < 2*half_size**2)
                & ( np.abs(galaxies_all["z_true"] - ztrue) < 0.017141 ) )          
            #vclu_gal = galaxies_all.loc[VICINITY_GAL]           # no need
        
            for clcl_gal in VICINITY_GAL:     
                vicinity_current_gal = galaxies_all.loc[clcl_gal]           
                vicenter_gal = list(zip(vicinity_current_gal["x_pix"].values*30-5, 
                                        vicinity_current_gal["y_pix"].values*30-5))
            
        if False: 
            sled = vicinity_current_gal[((vicinity_current_gal["x_pix"]*30-5 + 2.65)**2 + (vicinity_current_gal["y_pix"]*30-5 - 14.6)**2 < 0.1**2)]
            display(sled)
        
            print(sled["isub"].values[sled["isub"].values<100000])
            #plt.hist(sled["isub"].values[sled["isub"].values<100000])
            #plt.show()
            gall=[]
            for opopo in vicenter_gal:    
                #print(opopo)
                #print(np.sqrt((opopo[0] + 2.6)**2 + (opopo[1] - 14.5)**2))
                if ((opopo[0] + 2.65)**2 + (opopo[1] - 14.6)**2 < 0.1**2):
                    #print(opopo)
                    gall.append(opopo)
                    #display(vicinity_current_gal)            
            #print(clusters_all[clusters_all["ihal"] in sled["isub"].values]) 
            for s in sled["isub"].values:
                print("yes:"+str(s) if s in clusters_all["ihal"].values else "no")
                #print(s, clusters_all[clusters_all["ihal"]==s])
        
    if draw:
           
        trtr = plt.imshow(np.rot90(nmhg),
                          #convolve(np.rot90(nmhg), Gaussian2DKernel(2)), 
                          norm=matplotlib.colors.SymLogNorm(linthresh=LINTHRESH, linscale=1), 
                          origin='upper',
                          extent = axesforsmooth)     
        
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue', edgecolor='k', linewidth=1)
        plt.scatter(c_x_1, c_y_1, color='orangered', label = 'Centroid', edgecolor='k', linewidth=1)
        plt.scatter(xeeec, yeeec, color='dodgerblue', label = 'Max value', edgecolor='k', linewidth=1)
                       
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, 10*R, color='orangered', linestyle=":", lw=1, fill = False))
        
        if True:
            plt.gca().add_patch(plt.Circle(cntr, R*1.6, 
                               color='dodgerblue', linestyle="--", lw=2, fill = False))
            plt.gca().add_patch(plt.Circle(cntr, R*2.7, 
                               color='green', linestyle="--", lw=2, fill = False))
            plt.gca().add_patch(plt.Circle(cntr, R*8.1, 
                               color='grey', linestyle="--", lw=2, fill = False))
        
        x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
        y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
        y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]   
        plt.plot((x_s+5/60, x_s-5/60), (y_s, y_s), color='white')
        plt.text(x_s, y_S, f'10 arcmin $\\approx$ {10/60*D_A.value*np.pi/180:.0f} kpc', 
                 color='white', ha='center', va='center')
        
        plt.xlim(cntr[0]-half_size, cntr[0]+half_size)
        plt.ylim(cntr[1]-half_size, cntr[1]+half_size)
        plt.gca().set_aspect('equal', 'box')

        if current_cluster_number in clusters_34.index:
            plt.text(
                     (plt.gca().get_xlim()[1]-plt.gca().get_xlim()[0])*0.05+plt.gca().get_xlim()[0],
                     (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.05+plt.gca().get_ylim()[0],
                     f'M', color='white', ha='center', va='center', fontweight='bold')

#        if delete_superfluous:
            
#            for vv in vicenter:
#                plt.scatter(vv[0], vv[1], color='red', label = 'Subhaloes', s=3)
#                plt.gca().add_patch(plt.Circle((vv[0], vv[1]), vv[2], color='red', ls="-", lw=1, fill=False, alpha=0.5))        

           # for vv in vicenter_gal:
           #     plt.scatter(vv[0], vv[1], color='blue', label = 'Subhaloes', s=4)
            
            #for opopo in gall:
            #    plt.scatter(opopo[0], opopo[1], color='white', label = 'Subhaloes', s=3)
            #plt.gca().add_patch(plt.Circle((-2.65, 14.6), 0.1, color='white', ls="-", lw=1, fill=False))
                    
        plt.xlabel("RA, deg", size=13)
        plt.ylabel("DEC, deg", size=13)
        plt.xticks(size=13)
        plt.yticks(size=13)
        cb = plt.colorbar(trtr, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=13)
        
        if not ARF_weights:
            cb.set_label(f"Photons cm$^{{-2}}$ s$^{{-1}}$ arcmin$^{{-2}}$", size=13)
        else:
            cb.set_label(f"Counts s$^{{-1}}$ arcmin$^{{-2}}$", size=13)
            
        ttiittllee = f'#{current_cluster_number}: z={ztrue:.3f}, A={AREA:.1f} arcmin$^2, Res={ang_res:.1f}\'\'$'
        plt.gca().invert_xaxis()
        
        #print(plt.gca().get_xlim()*u.deg)
        #print(plt.gca().get_ylim()*u.deg)
        
        plt.title(ttiittllee, fontsize=13)
        #cb.ax.set_yticklabels(['0', '1', '10', '100', '1000'])
        
        handles, labels = plt.gca().get_legend_handles_labels()
        #l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label=str(r)+"$\\cdot R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l2])
        #plt.legend(handles=handles, loc=3, fontsize=13)
        #plt.show()

    if delete_superfluous:
        return np.fliplr(np.rot90(nmhg_unfiltered)), np.fliplr(np.rot90(1-nmhg_mask))
    else:
        return np.fliplr(np.rot90(nmhg)) #*(wedge(1)+wedge(2)+wedge(3)+wedge(4))


def create_circle_mask(x_center, y_center, radius, N):
    # Initialize an NxN mask with zeros
    mask = np.zeros((N, N), dtype=np.uint8)

    # Create a grid of coordinates
    Y, X = np.ogrid[0:N, 0:N]

    # Calculate the distance from the center of the circle
    dist_from_center = np.sqrt((X - x_center)**2 + (Y - y_center)**2)

    # Create a mask for the circle
    mask[dist_from_center <= radius] = True  # Set pixels within the circle to True
    
    return mask


def koltso(r_in, r_out, mm, NMHG, d_pixels):
    # Create a full-sized mask with the same shape as NMHG
    mask = np.zeros(NMHG.shape, dtype=bool)
    
    # Define the limits for the slice
    x_min = max(mm[0] - d_pixels, 0)
    x_max = min(mm[0] + d_pixels + 1, NMHG.shape[0])
    y_min = max(mm[1] - d_pixels, 0)
    y_max = min(mm[1] + d_pixels + 1, NMHG.shape[1])
    
    # Extract the relevant slice from NMHG
    kusok = NMHG[x_min:x_max, y_min:y_max]
    
    #plt.imshow(kusok)
    #plt.show()
    
    # Create a grid of distances from the center
    Y, X = np.indices((x_max - x_min, y_max - y_min))
    dist_from_center = np.sqrt((X - d_pixels)**2 + (Y - d_pixels)**2)
    
    # Create the mask for the circular area
    circular_mask = (dist_from_center <= r_out) & (dist_from_center >= r_in)

    # Place the circular mask back into the full-sized mask
    mask[x_min:x_max, y_min:y_max] = circular_mask
    
    # Return the masked values and the mask
    return NMHG * mask, mask
    
    
def wedge(n, l=201):

    w = np.zeros((l,l))
    
    if n == 1:
        w[:l//2, l//2:] = 1
    elif n == 2:
        w[:l//2, :l//2] = 1
    elif n == 3:
        w[l//2:, :l//2] = 1
    elif n == 4:
        w[l//2:, l//2:] = 1
    else:
        print("?")
     
    return w


    
def brightness_profile(clusternumber, hist, mmmask, field_length, draw=True, ARF_weights=False, errors=False):
    
#    print(clusternumber)
    cc = clusters.loc[clusternumber]
    
    snap_id_str = binned_clusters[clusternumber][1]   # id of photon list

    if snap_id_str == '124':
        zslice = 0.174192889973847
    elif snap_id_str == '128':
        zslice = 0.13708140389145
    elif snap_id_str == '132':
        zslice = 0.101142861718869
    elif snap_id_str == '136':
        zslice = 0.0663401914452304
    elif snap_id_str == '140':
        zslice = 0.032637492755919
    
    R_500 = cc["R500"]/0.704/(1+zslice)  # kpc
    ztrue = cc["z_true"]
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    R_500_rescaled = R_500/D_A.value*180/np.pi*60      # arcmin
#    print(R_500_rescaled)
    #ang_res = 2*10*R_500_rescaled*60 / len(hist)    # arseconds
#    print(ang_res)
    
    r_pixels_max = int(len(hist)/2)                  # depends on field size
    r500r = int(r_pixels_max/(field_length/2))       # field length should be in units of R500
    setka_bins = np.append([0, 2, 3, 4], 
                           np.geomspace(5, r_pixels_max, 50)) # .astype(int)       # borders of bins
    setka = [(a+b)/2 for a, b in zip(setka_bins[:-1], setka_bins[1:])]             # centers of bins
    c2 = [r_pixels_max, r_pixels_max]                # center of field
    err = np.diff(setka_bins)/2                      # just bins width
    brightness = []
    br1, br2, br3, br4 = [], [], [], []              # in 4 wedges    
    
    for i in tqdm(range(0, len(setka))):
        
        ring = koltso(setka_bins[i], setka_bins[i+1], c2, hist, r_pixels_max-1)       # returns both ring and mask for it
        
        # test place
        
        if True and (i>30) and (i<45):
            
            plt.figure(figsize=(22, 6))
            
            aop = sum(ring[1].flatten())
            print(aop)
            pixsize = 1       # R500 = 10 pix = 10 arcmin => 1 arcmin = 1 pix
            
            ring0 = ring[0]
            ring0 = ring0 * 10000 * aop
                        
            plt.subplot(131)
            plt.imshow(ring0, norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.imshow(hist, norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1), alpha=0.5)
            half_length = int(len(reduced_hist)/2)
            r500r = int(len(reduced_hist)/20)
            plt.gca().add_patch(plt.Circle((half_length, half_length), r500r, 
                                color='orangered', linestyle="--", lw=2, fill = False))
            plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*1.6, 
                                color='dodgerblue', linestyle="--", lw=2, fill = False))
            plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*2.7, 
                                color='green', linestyle="--", lw=2, fill = False))
            plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*8.1, 
                                color='grey', linestyle="--", lw=2, fill = False))
            cb = plt.colorbar(fraction=0.046, pad=0.04)
            cb.set_label(f"Counts s$^{{-1}}$ arcmin$^{{-2}}$", size=13)                       
            skolko = 11
        #    plt.xticks(np.linspace(0, 1, skolko)*20, np.linspace(-10, 10, skolko).astype(int))
        #    plt.yticks(np.linspace(0, 1, skolko)*20, np.linspace(10, -10, skolko).astype(int))
            plt.axvline(int(len(reduced_hist)/2), color='white', ls='--', lw=2)
            plt.axhline(int(len(reduced_hist)/2), color='white', ls='--', lw=2)            
            
            plt.subplot(132)
            rrrr = [(ring0)[ii, jj] for ii in range(len(reduced_hist)) for jj in range(len(reduced_hist)) if (ring[1])[ii, jj] == 1]
            rrrr = np.array(rrrr)
            #rrrr = rrrr[rrrr>0]
            rrrr[rrrr==0] = 1e-7
            #print(rrrr)
            srr = np.mean(rrrr)
            RMS = np.std(rrrr) # np.sqrt(sum([(el**2) for el in rrrr])/len(rrrr))
            VAR = np.var(rrrr)
            beans = np.logspace(np.log10(min(rrrr)), np.log10(max(rrrr)), 20)
            beans = np.linspace(0, 500, 20)           
            Vals, Bins, _ = plt.hist(rrrr, bins=beans, alpha=0.7) 
                            #, density=False, label=f'mean = {srr:.7f}\nrms = {RMS:7f}')
            BinsC = [(a+b)/2 for a, b in zip(Bins[:-1], Bins[1:])]
            BinsD = np.diff(Bins)
            plt.xlabel("Counts", size=13) #  s$^{{-1}}$ arcmin$^{{-2}}$
            #plt.scatter(BinsC, Vals, label=f'mean = {srr:.7f}\nvar = {VAR:5f}') # if np.histogram
            plt.axvline(srr, color='red', ls='--', lw=3, label=f'mean = {srr:.2f}')# $\\times 10^{{-5}}$')
            #plt.axvline(srr-RMS, color='green', ls='--', lw=2)
            #plt.axvline(srr+RMS, color='green', ls='--', lw=2, label=f'rms = {RMS:7f}')
     #       plt.axvline(np.dot(BinsC,Vals)/sum(Vals), 
     #                   color='red', ls='--', lw=5, alpha=0.5, label='sum(bins_centers*vals)/sum(vals)')
           # plt.xscale("log")
            #plt.yscale("log")
            
      #      print(srr)
      #      ps1 = stats.poisson.rvs(mu=srr/20, size=aop)*12
      #      plt.hist(ps1, bins=beans, weights=np.ones_like(ps1, dtype=float), alpha=0.5)
      #      ps2 = stats.poisson.rvs(mu=srr/20, size=aop)*15
      #      vvv, bbb, _ = plt.hist(ps2+ps1, bins=beans, weights=np.ones_like(ps2, dtype=float), alpha=0.4)
      #      plt.plot(bbb[:-1], vvv)
            
            n_processes = 27  # 12+15=27 unit processes
            ps_unit = np.sum([stats.poisson.rvs(mu=srr/20/27, size=aop) for _ in range(n_processes)], axis=0)
            ps_unit = stats.poisson.rvs(mu=srr/20, loc=0, size=aop)
            plt.hist(ps_unit*25, bins=np.linspace(0, 500, 20), density=False, weights=np.ones_like(ps_unit, dtype=float), alpha=0.5)
            
        #    r = rrrr.astype(int)
        #    print(r)
        #    print(r.mean())
        #    k_min = int(r.min())
        #    k_max = int(r.max())
        #    bins = np.arange(k_min - 0.5, k_max + 1.5, 1)
        #    plt.hist(r, bins=bins, density=True, alpha=0.4, label='data')
        #    k = np.arange(k_min, k_max + 1)
        #    pmf = stats.poisson(mu=r.mean()/2).pmf(k)
        #    plt.plot(k, pmf*200, 'o-', label=fr'Poisson($\lambda$={r.mean():.2f})')
            
                        
            wedges = [0, 0, 0, 0]
            for we in [1,2,3,4]:
                cw1      = (ring[0])*wedge(we)
                cw1_ring = (ring[1])*wedge(we)
                ww = [cw1[ii, jj] for ii in range(len(reduced_hist)) for jj in range(len(reduced_hist)) if cw1_ring[ii, jj] == 1]
                ww = np.array(ww)
                ww[ww==0] = 1e-7
                mw = np.mean(ww)
                wedges[we-1] = mw            
            
        #    print(wedges)
            std_we = np.std(wedges)*10000*aop
        #    print(std_we)
            
            plt.axvline(srr+std_we, color='orange', ls='--', lw=3, label=f'std (4)  = {std_we:.2f}')# $\\times 10^{{-5}}$')           
            plt.axvline(srr-std_we, color='orange', ls='--', lw=3, label=f'mean $\\pm$ std')
            plt.axvline(srr+RMS, color='magenta', ls='--', lw=3, label=f'std all  = {RMS:.2f}')# $\\times 10^{{-5}}$')           
            plt.axvline(srr-RMS, color='magenta', ls='--', lw=3, label=f'mean $\\pm$ std')
            plt.legend()
   #         plt.xlim(1e-7, 1e-4)
            
            plt.subplot(133)
    #        print(np.sum(np.array(Vals[:])*np.array(BinsC[:])), sum(ring[0].flatten()))
    #        print((hist*ring[1]).sum()/sum(ring[1].flatten()))
            for v in range(0, len(Vals)):
                #yy = np.dot(np.array(Vals[:-v]), np.array(BinsC[:-v])
                filterrrred = rrrr[rrrr <= BinsC[-v]]
                yy = np.mean(filterrrred)
                plt.scatter(BinsC[-v], yy, color='k')
            plt.axhline(srr, color='red', ls='--', lw=2)
            plt.xlabel("Counts s$^{{-1}}$ arcmin$^{{-2}}$ (right end cutoff)", size=13)
            plt.ylabel("Mean value of surface brightness", size=13)       
            #x1 = np.arange(0, max(rrrr)+1)
            #st = stats.poisson.pmf(x1, int(srr))
            #plt.plot(x1, st/max(st), 'o-', label=f'Poisson(λ={srr:.5f})')
            plt.xscale("log")
            plt.yscale("log")
            #plt.legend()
            #plt.tight_layout()
   #         plt.savefig(str(i)+'.png', dpi=200)
            plt.show()
            plt.figure(figsize=(6, 6))
                
        ####     hist*ring[1] = ring[0] ---- fact
        
        if isinstance(mmmask, str):
            nbvc = np.array([0.])
            mmmask = np.ones_like(hist)
        else:
            nbvc = ring[1]*(1-mmmask)            # reduce the area of mask for ring - only for individual profiles!
        
        cheese = ring[1] - nbvc    # actually ring[1]*mmask        
        
        if False:
            print('Total brightness inside the ring:', ring[0].sum())
            print('Total area of ring in pixels:', sum(ring[1].flatten()))
            print('Area to exclude:', sum(nbvc.flatten()))
            print('Total brightness after filtering:', sum( (hist*cheese).flatten() ) )
            print('Total area of ring excluding masked regions:', sum(ring[1].flatten())-sum(nbvc.flatten()))
            
            plt.figure(figsize=(8,8))
            plt.subplot(221)
            plt.imshow(ring[1]+3*(1-mmmask), origin='upper')# norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(222)
            plt.imshow(ring[0], origin='upper', norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(223)
            plt.imshow(cheese, origin='upper')#, norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(224)
            plt.imshow(hist*cheese, origin='upper', norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.show()
            plt.figure(figsize=(6, 6))
        
        #print(ring[0].sum(), sum(ring[1].flatten()))
        #print()
        
        br = (hist*cheese).sum()/sum(cheese.flatten())
        
        #ring[0].sum()/(sum(ring[1].flatten())-sum(nbvc.flatten()))      # (pix_in_k2-pix_in_k1)
        
        brightness.append(br)
             
        if errors:
            cw = cheese*wedge(1)
            brw1 = (hist*cw).sum()/sum(cw.flatten())
            cw = cheese*wedge(2)
            #print(cw)
            brw2 = (hist*cw).sum()/sum(cw.flatten())
            #print(brw2)
            cw = cheese*wedge(3)
            brw3 = (hist*cw).sum()/sum(cw.flatten())
            cw = cheese*wedge(4)
            brw4 = (hist*cw).sum()/sum(cw.flatten())
            
            br1.append(brw1)
            br2.append(brw2)
            br3.append(brw3)
            br4.append(brw4)
            
    #print(br1, br2, br3, br4)
    
    meanbr = [np.mean([a, b, c, d]) for a,b,c,d in zip(br1, br2, br3, br4)]
    stdbr = [np.std([a, b, c, d]) for a,b,c,d in zip(br1, br2, br3, br4)]
    
    #print(meanbr, stdbr)
    #print(brightness)
    #print(np.array(setka)/r500r*(10*998/1000))
    #print(len(setka), len(brightness))
    
    # Actually R500inmin should be since 10*998/1000  
    # this is because of 343 Mpc: 10 / (1/343*180/pi*60) = 10 / 10.022 = 0.998 
    
    R500inmin = 10 #R_500_rescaled

    rr = np.array(setka)/r500r*R500inmin 
    dr = np.diff(setka_bins)/r500r*R500inmin
    
    #print(rr, dr)
    
    if draw:
        
        plt.xlabel("Radius, arcmin", fontsize=12)  # "Radius in units of $R_{500}$")
        
        if not ARF_weights:
            plt.ylabel("Photons cm$^{{-2}}$ s$^{{-1}}$ arcmin$^{{-2}}$", fontsize=12) # "Brightness in relative units")
        else:
            plt.ylabel("Counts s$^{{-1}}$ arcmin$^{{-2}}$", fontsize=12) # "Brightness in relative units")
        
        plt.xscale("log")
        plt.yscale("log")

        plt.axvline(R500inmin, linestyle='--', color='orangered')#, label='$R_{500c}$', lw=1)
        plt.axvline(R500inmin*1.6, linestyle='--', color='dodgerblue')#, label='$R_{200c} = 1.6 \\cdot R_{500c}$', lw=1)
        plt.axvline(R500inmin*2.7, linestyle='--', color='green')#, label='$R_{200m} = 2.7 \\cdot R_{500c}$', lw=1)
        plt.axvline(R500inmin*8.1, linestyle='--', color='grey')#, label='$R_{ta} = 8.1 \\cdot R_{500c}$', lw=1)
        
#        plt.text(19, 1e-1, '$R_{500c}:R_{200c}:R_{200m}:R_{ta}=$\n$=1:1.6:2.7:8.1$', fontsize=11,
#                 bbox=dict(facecolor='white', alpha=0.99, edgecolor='grey'),
#                 ha='center', va='center')        
        #plt.scatter(setka[:-1]/r500r*(10*998/1000), np.array(brightness)/10000, color='black', s=7)
        #print(np.array(setka)/r500r*R500inmin)
                
        plt.errorbar(rr, #                       np.array(setka)/r500r*R500inmin, 
                     np.array(brightness), 
                     xerr=err/r500r*R500inmin, linewidth=0, marker='o', markersize=3, alpha=0.95,
                     elinewidth=1, capsize=0, color='black', label='Stacked image')
        #plt.ylim(1e-5, 5e-1)
        plt.legend(loc=3, fontsize=12)
        plt.xticks([0.1, 1, 10, 100], [0.1, 1, 10, 100])
        #plt.gca().set_aspect('auto', 'box')
        #plt.show()
    
        if errors:
            plt.errorbar(rr, meanbr, yerr=stdbr, lw=0,
                         fmt='', capsize=0, capthick=1, elinewidth=1, color='black', ecolor='black', alpha=0.95)
        
        # 4 different plot for 4 wedges:
        if False:                 
            plt.plot(rr, np.array(br1))
            plt.plot(rr, np.array(br2))
            plt.plot(rr, np.array(br3))
            plt.plot(rr, np.array(br4))
        
        resc1 = lambda x: x/R500inmin
        resc2 = lambda x: x*R500inmin
        
        ax2 = plt.gca().secondary_xaxis("top", functions=(resc1, resc2))
        ax2.set_xlabel("Radius / R$_{500}$", fontsize=12)
        ax2.set_xticks([0.01, 0.1, 1, 10], [0.01, 0.1, 1, 10])
        
    # calculating luminosity
    
    D_L = FlatLambdaCDM(H0=100*0.704, Om0=0.272).luminosity_distance(ztrue).value # Mpc
    
    N = np.sum(2*np.pi*rr[:16]*dr[:16]*brightness[:16])  # 16th position corresponds to ~R500
    
    #N = np.sum(2*np.pi*rr*dr*brbr)
    #print(N)
    flx = N * (1.1*1.60218e-9) / 140
    #print(flx)
    L_p = flx * 4 * np.pi * D_L**2 * (3.08e24)**2
    #print(L_p)
    #L_catalogue = clusters.loc[clusternumber]["Lx500"],'e+44'
    #print(L_spec)
   
    return np.array(setka)/r500r*R500inmin, np.array(brightness), R_500_rescaled/10, L_p
    

def draw_84_panels():

    NNN = 84
    
    size = 6

    plt.figure(figsize=((size)*7+6*3, 5*12+11*2.5))
    plt.tight_layout()
    
    for cl_num in tqdm(clusters.index[:NNN]):
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        pho_hist = extract_photons_from_cluster(cl_num, draw=True, delete_superfluous=True, histlen=201)


def calc_l_T(T, T_left, T_right, Xplot=False):
  
    x.Xset.chatter = 0
    x.AllData.clear()
    x.AllData.removeDummyrsp()
    x.AllData.dummyrsp(lowE=0.1, highE=10.0, nBins=1024)
    x.Xset.addModelString("APEC_TRACE_ABUND", "0")
    x.Xset.abund = "angr"

    if Xplot:
        x.Plot.device = '/xs'
    else:
        x.Plot.device = '/null'
    
    expo = 40000
    Ab = 0.2
    z = 0.0
    Norm = 1/4/np.pi/(1+z)**2
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
    ARF_NAME = '../erosita/tm1_arf_filter_000101v02.fits'
    
    #RMF_NAME = 'rmf_v1.fits'
    #ARF_NAME = 'esf10.Dsur1234regR3cCaXv2.0001.arf'
        
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
    
    
def E(z):
    
    O_M, O_L, O_K = 0.272, 0.728, 0.000
    
    return np.sqrt( O_M*(1+z)**3 + O_K*(1+z)**2 + O_L )
    

def koltso_old(r_in, r_out, mm, NMHG, d_pixels):

    kusok = np.zeros((2*d_pixels+1, 2*d_pixels+1))
        
    for i in range(mm[0]-d_pixels, mm[0]+d_pixels+1):
        for j in range(mm[1]-d_pixels, mm[1]+d_pixels+1):
            kusok[i - (mm[0]-d_pixels)][j - (mm[1]-d_pixels)] = NMHG[i][j]
        
    Y, X = np.ogrid[(mm[1]-d_pixels):(mm[1]+d_pixels+1), (mm[0]-d_pixels):(mm[0]+d_pixels+1)]
    dist_from_center = np.sqrt((X - mm[0])**2 + (Y-mm[1])**2)
    
    mask = (dist_from_center <= r_out) & (dist_from_center >= r_in)
    
    #print(kusok[mask])
    #plt.subplot(121)
    #plt.imshow(mask[990:1010, 990:1010])
    #plt.subplot(122)
    #plt.imshow((kusok*mask)[990:1010, 990:1010])
    #plt.show()
        
    return mask*kusok, mask
 

def kruzhok_old(r_pixels, mm, NMHG, d_pixels):

    kusok = np.zeros((2*d_pixels+1, 2*d_pixels+1))
        
    for i in range(mm[0]-d_pixels, mm[0]+d_pixels+1):
        for j in range(mm[1]-d_pixels, mm[1]+d_pixels+1):
            kusok[i - (mm[0]-d_pixels)][j - (mm[1]-d_pixels)] = NMHG[i][j]
        
    Y, X = np.ogrid[(mm[1]-d_pixels):(mm[1]+d_pixels+1), (mm[0]-d_pixels):(mm[0]+d_pixels+1)]
    dist_from_center = np.sqrt((X - mm[0])**2 + (Y-mm[1])**2)
    
    mask = dist_from_center <= r_pixels
        
    return mask*kusok, mask
