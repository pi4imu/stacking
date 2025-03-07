def draw_all_profiles():

    r_pixels_max = 1000 # 5*r500r    # depends on field size
    r500r = int(r_pixels_max/10)
    setka_bins = np.append([0, 1, 2, 3, 4],np.geomspace(5, r_pixels_max, 20))
    setka = [(a+b)/2 for a, b in zip(setka_bins[:-1], setka_bins[1:])]  # centers of bins
    err = np.diff(setka_bins)/2

    yar = np.array(yarkosti.to_numpy()[8])*0

    for i in range(0, len(yarkosti[::])):
    
        one = yarkosti.to_numpy()[i]
        
        if one[0]==0 or one[1]==0:
            #print(one)
            llww=1
        else:
            llww=1
    
        f1 = 1000/R500S[i]
        f2 = E(reds[i])**(-4)*(1+reds[i])**3
        factor = f1*f2
        one = one*factor
    
        yar = yar + one
    
        plt.plot(np.array(setka)/r500r*(10*998/1000), 
                 np.array(one), 
                 linewidth=llww, marker='.', markersize=0, alpha=0.75,
                 color='black')
    
    plt.plot(np.array(setka)/r500r*(10*998/1000), 
             np.array(yar)/84, 
             linewidth=3, marker='.', 
             markersize=0, alpha=0.95,
             color='orangered')

    #plt.errorbar(np.array(setka)/r500r*(10*998/1000), 
    #             np.array(yar)/84, 
    #             xerr=err/r500r*(10*998/1000), 
    #             linewidth=0, marker='.', 
    #             markersize=3, alpha=0.95,
    #             elinewidth=1, capsize=0, color='orangered')

    plt.xlabel("Radius, arcmin", fontsize=12)  # "Radius in units of $R_{500}$")
    plt.ylabel("Counts s$^{{-1}}$ arcmin$^{{-2}}$", fontsize=12) # "Brightness in relative units")

    plt.xscale("log")
    plt.yscale("log")

    plt.axvline(10*998/1000, linestyle='--', color='orangered', label='$R_{500c}$', lw=2)
    plt.axvline(10*998/1000*1.6, linestyle='--', color='dodgerblue', label='$R_{200c} = 1.6 \\cdot R_{500c}$', lw=2)
    plt.axvline(10*998/1000*2.7, linestyle='--', color='green', label='$R_{200m} = 2.7 \\cdot R_{500c}$', lw=2)
    plt.axvline(10*998/1000*8.1, linestyle='--', color='magenta', label='$R_{ta} = 8.1 \\cdot R_{500c}$', lw=2)

    #plt.ylim(1e-8, 4e-2)

    plt.legend(loc=3, fontsize=12)
    plt.xticks([0.1, 1, 10, 100], [0.1, 1, 10, 100])
    plt.gca().set_aspect('auto', 'box')
       
    return None
