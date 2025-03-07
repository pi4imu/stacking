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
    
    # D_A = 343000
    
    #R_500_fid = 1000/343000*180/np.pi   # degrees
            
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
    
   # t1 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_124.fits", hdu=2)
   # t2 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_128.fits", hdu=2)
   # t3 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_132.fits", hdu=2)
   # t4 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_136.fits", hdu=2)
   # t5 = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_140.fits", hdu=2)    
   # t = vstack([t1, t2, t3, t4, t5])
    
    if withagn:
        for snapid in [124, 128, 132, 136, 140]:
            t1 = Table.read("../data/eROSITA_30.0x30.0/Phox/AGNphlist_"+str(snapid)+".fits", hdu=2)
            t = vstack([t, t1])
    
    SLICE = t.to_pandas()        # for photons extraction
    SLICE1 = t.to_pandas()       # for drawing
    SLICE2 = t.to_pandas()       # for center searching if it is not from table
    #SLICE3 = t.to_pandas()       # for rescaling SLICE2
    
    R = r * R_500_rescaled
        
    AREA = np.pi*R**2*3600   # min2
    
    #setting area and resolution for searching for center
    
    ang_res = 4                    # this value will be changed later 
    halfsidelength = 10           # in units of R500
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
            hs4s = half_size/3/(halfsidelength/3)
        else:
            hs4s = half_size/1/(halfsidelength/3)
            
        # making 2D histogram with side length 2*half_size with center (RA_c, DEC_c) without drawing
                
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < hs4s) & (np.abs(SLICE2["DEC"]-DEC_c) < hs4s), True, False)
        search4centroid = SLICE2[SLICE2['what'] == True]
        search4centroid = search4centroid.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(search4centroid["RA"], search4centroid["DEC"], bins=int(2*hs4s*3600/ang_res))
                       
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
        
        SLICE["check"]=np.where((SLICE["RA"]-c_x_1)**2+(SLICE["DEC"]-c_y_1)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
    
                   
    # this goes to final panels image centered at cntr (which is set to centroid position by default)
           
    SLICE1["todraw"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size),
                                    True, False)
    SLICE1 = SLICE1[SLICE1['todraw'] == True]
    SLICE1 = SLICE1.drop("todraw", axis=1)
        
    SLICE1 = SLICE1[SLICE1['ENERGY']>0.3]           # previosly 0.3-2.3
    SLICE1 = SLICE1[SLICE1['ENERGY']<2.3]
    
    
    # drawing without or with applying weights that are given by ARF
    
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
    #factor=1
    #print(factor)

    nmhg = nmhg*factor
        
    # deleting less massive haloes
        
    if delete_superfluous:
        
        # taking haloes only from the slice to which this cluster belongs 
        
        VICINITY = np.where( ((clusters_all["x_pix"]*30-5 - c_x_1)**2 + (clusters_all["y_pix"]*30-5 - c_y_1)**2 < 2*half_size**2) & ( np.abs(clusters_all["z_true"] - ztrue) < 0.017141 ) )
        #print(VICINITY)
        vclu = clusters_all.loc[VICINITY]  
        #display(vclu)
                
        for clcl in VICINITY:
        
            vicinity_current = clusters_all.loc[clcl].drop(current_cluster_number)

            D_A_ = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(vicinity_current["z_true"])*1000 # kpc
            R_500_ = vicinity_current["R500"]*0.704
            R_500_s = R_500_/D_A_*180/np.pi # degrees
            
            vicenter = list(zip(vicinity_current["x_pix"].values*30-5, 
                                vicinity_current["y_pix"].values*30-5, 
                                R_500_s))
            #print(clcl)
            #print(vicenter)
            
    # postponed attempt to take galaxies into account        
            
    #    galaxies_all = pd.read_csv("../data/eROSITA_30.0x30.0/Catalouges/galaxies.dat", sep='\s+', header=0)    
    #    VICINITY_GAL = np.where( ((galaxies_all["x_pix"]*30-5 - c_x_1)**2 + (galaxies_all["y_pix"]*30-5 - c_y_1)**2 < 2*half_size**2) & ( np.abs(galaxies_all["z_true"] - ztrue) < 0.017141 ) )          
    #    vclu_gal = galaxies_all.loc[VICINITY_GAL]          
    #    for clcl_gal in VICINITY_GAL:     
    #        vicinity_current_gal = galaxies_all.loc[clcl_gal]           
    #        vicenter_gal = list(zip(vicinity_current_gal["x_pix"].values*30-5, vicinity_current_gal["y_pix"].values*30-5))
    
        # creating the mask which will go to output along with unfiltered field
        
        nmhg_mask = np.zeros((histlen, histlen))
        MASK = pd.DataFrame([])
        
        for vic in vicenter:
        
            #print(vic)
            ra_pix = np.abs(vic[0]-cntr[0]-half_size)*3600/ang_res
            dec_pix = (vic[1]-cntr[1]+half_size)*3600/ang_res
            radius_pix = vic[2]*3600/ang_res
            
            ra_pix = ra_pix.astype(int)
            dec_pix = dec_pix.astype(int)
            radius_pix = radius_pix.astype(int)
            
            circle_mask = create_circle_mask(ra_pix, dec_pix, radius_pix, len(nmhg_mask))
            if circle_mask[int(histlen/2), int(histlen/2)] == 0:
                nmhg_mask = nmhg_mask + circle_mask
        
        if False:
            pup = create_circle_mask(400, 700, 333, 2001)
            nmhg_mask = nmhg_mask + pup
                
        nmhg_mask[nmhg_mask > 1] = True   
        nmhg_mask = np.rot90(nmhg_mask)         # some magic
        
        nmhg_unfiltered = nmhg
        nmhg = nmhg_unfiltered * (1-nmhg_mask)
                
    
    if draw:
           
        trtr = plt.imshow(np.rot90(nmhg),
                          #convolve(np.rot90(nmhg), Gaussian2DKernel(2)), 
                          norm=matplotlib.colors.SymLogNorm(linthresh=LINTHRESH, linscale=1), 
                          origin='upper',
                          extent = axesforsmooth)     
        
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue')
        plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
                        
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, 10*R, color='orangered', linestyle=":", lw=1, fill = False))
        
        x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
        y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
        y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]   
        plt.plot((x_s+5/60, x_s-5/60), (y_s, y_s), color='white')
        plt.text(x_s, y_S, f'10 arcmin $\\approx$ {10/60*D_A.value*np.pi/180:.0f} kpc', 
                 color='white', ha='center', va='center')
        
        plt.xlim(cntr[0]-half_size, cntr[0]+half_size)
        plt.ylim(cntr[1]-half_size, cntr[1]+half_size)
        plt.gca().set_aspect('equal', 'box')
        
        if delete_superfluous:
                 
         #   for vv in vicenter_gal:
         #       plt.scatter(vv[0], vv[1], color='blue', label = 'Subhaloes', s=7)
                                                       
            for vv in vicenter:
                plt.scatter(vv[0], vv[1], color='red', label = 'Subhaloes', s=3)
                plt.gca().add_patch(plt.Circle((vv[0], vv[1]), vv[2], color='red', ls="-", lw=1, fill=False))        
        
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
        return nmhg_unfiltered, 1-nmhg_mask
    else:
        return nmhg #*(wedge(1)+wedge(2)+wedge(3)+wedge(4))


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
    
    
def wedge(n, l=2001):

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
                           np.geomspace(5, r_pixels_max, 20)) # .astype(int)       # borders of bins
    setka = [(a+b)/2 for a, b in zip(setka_bins[:-1], setka_bins[1:])]             # centers of bins
    c2 = [r_pixels_max, r_pixels_max]                # center of field
    err = np.diff(setka_bins)/2                      # just bins width
    brightness = []
    br1, br2, br3, br4 = [], [], [], []              # in 4 wedges    
    
    for i in tqdm(range(0, len(setka))):
        
        ring = koltso(setka_bins[i], setka_bins[i+1], c2, hist, r_pixels_max-1)       # returns both ring and mask for it
        
        ####     hist*ring[1] = ring[0] ---- fact
        
        if isinstance(mmmask, str):
            nbvc = np.array([0.])
            mmmask = np.ones_like(hist)
        else:
            nbvc = ring[1]*(1-mmmask)            # reduce the area of mask for ring 
        
        cheese = ring[1] - nbvc        
        
        if False:
            print('Total brightness inside the ring:', ring[0].sum())
            print('Total area of ring in pixels:', sum(ring[1].flatten()))
            print('Area to exclude:', sum(nbvc.flatten()))
            print('Total brightness after filtering:', sum( (hist*cheese).flatten() ) )
            print('Total area of ring excluding masked regions:', sum(ring[1].flatten())-sum(nbvc.flatten()))
            
            plt.figure(figsize=(8,8))
            plt.subplot(221)
            plt.imshow(ring[1]+3*(1-mmmask), origin='lower')# norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(222)
            plt.imshow(ring[0], origin='lower', norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(223)
            plt.imshow(cheese, origin='lower')#, norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
            plt.subplot(224)
            plt.imshow(hist*cheese, origin='lower', norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1))
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

        plt.axvline(R500inmin, linestyle='--', color='orangered', label='$R_{500c}$', lw=1)
        plt.axvline(R500inmin*1.6, linestyle='--', color='dodgerblue', label='$R_{200c} = 1.6 \\cdot R_{500c}$', lw=1)
        plt.axvline(R500inmin*2.7, linestyle='--', color='green', label='$R_{200m} = 2.7 \\cdot R_{500c}$', lw=1)
        plt.axvline(R500inmin*8.1, linestyle='--', color='magenta', label='$R_{ta} = 8.1 \\cdot R_{500c}$', lw=1)

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
            plt.errorbar(rr, meanbr, yerr=stdbr,
                         fmt='.', capsize=0, capthick=1, elinewidth=1, color='black', ecolor='black')
        
        # 4 different plot for 4 wedges:
        if True:                 
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
   
    return np.array(setka)/r500r*R_500_rescaled, np.array(brightness), R_500_rescaled/10, L_p
    

def draw_84_panels():

    NNN = 84
    
    size = 6

    plt.figure(figsize=((size)*7+6*3, 5*12+11*2.5))
    plt.tight_layout()
    
    for cl_num in tqdm(clusters.index[:NNN]):
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        pho_hist = extract_photons_from_cluster(cl_num, draw=True, delete_superfluous=False)


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
    
    
def draw_total_hist(hhhh):
  
    length = len(hhhh)

    half_length = int(length/2)
    r500r = int(half_length/10)

    plt.imshow(np.rot90(hhhh), norm=matplotlib.colors.SymLogNorm(linthresh=0.0001, linscale=1), origin='upper')
    plt.colorbar(fraction=0.046, pad=0.04)
    #cb.set_label(f"Counts s$^{{-1}}$ arcmin$^{{-2}}$", size=13)

    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r, 
                                color='orangered', linestyle="--", lw=3, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*1.6, 
                                color='dodgerblue', linestyle="--", lw=3, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*2.7, 
                                color='white', linestyle="--", lw=3, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*8.1, 
                                color='grey', linestyle="--", lw=3, fill = False))

    x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
    y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
    y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]   
    plt.plot((x_s+r500r/2, x_s-r500r/2), (y_s, y_s), color='white')
    plt.text(x_s, y_S, f'10 arcmin $\\approx$ 1 Mpc', 
            color='white', ha='center', va='center')

    plt.xlabel("$20 \\times R_{500}$", fontsize=12)
    plt.ylabel("$20 \\times R_{500}$", fontsize=12)

    plt.show()
