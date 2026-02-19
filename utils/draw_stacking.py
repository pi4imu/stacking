def draw_stacked_image(histogram, r500r):

    plt.imshow(histogram, norm=matplotlib.colors.SymLogNorm(linthresh=0.000001, linscale=1), origin='upper')
    
    cb = plt.colorbar(fraction=0.046, pad=0.04)
    cb.set_label(f"Counts s$^{{-1}}$ arcmin$^{{-2}}$", size=13)

    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r, 
                               color='red', linestyle="--", lw=2, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*1.6, 
                               color='cyan', linestyle="--", lw=2, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*2.7, 
                               color='mediumspringgreen', linestyle="--", lw=2, fill = False))
    plt.gca().add_patch(plt.Circle((half_length, half_length), r500r*8.1, 
                               color='gainsboro', linestyle="--", lw=2, fill = False))

    x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
    y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.92+plt.gca().get_ylim()[0]
    y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]   
    plt.plot((x_s+r500r/2, x_s-r500r/2), (y_s, y_s), color='white', lw=2)
    plt.text(x_s, y_S, f'10 arcmin $\\approx$ 1 Mpc', color='white', ha='center', va='center',
             #bbox=dict(facecolor='black', alpha=0.1), 
             fontweight='extra bold')
    
    skolko = 11
    plt.xticks(np.linspace(0, 1, skolko)*20*r500r, np.linspace(-10, 10, skolko).astype(int))
    plt.yticks(np.linspace(0, 1, skolko)*20*r500r, np.linspace(10, -10, skolko).astype(int))

    plt.xlabel("$ (x-x_c) \\ / \\ R_{500}$", fontsize=13)
    plt.ylabel("$ (y-y_c) \\ / \\ R_{500}$", fontsize=13)
    
#    plt.gca().add_patch(plt.Circle((half_length, half_length), 10*r500r, color='orangered', linestyle=":", lw=1, fill = False))
    
    #plt.axvline(0, linestyle='--', color='orangered', label='$R_{500c}$', lw=1)
    #plt.axvline(0, linestyle='--', color='dodgerblue', label='$R_{200c} = 1.6 \\cdot R_{500c}$', lw=1)
    #plt.axvline(0, linestyle='--', color='green', label='$R_{200m} = 2.7 \\cdot R_{500c}$', lw=1)
    #plt.axvline(0, linestyle='--', color='grey', label='$R_{ta} = 8.1 \\cdot R_{500c}$', lw=1)
    #plt.legend(fontsize=11, frameon=True, loc=3)
        
    return None
    
def draw_stacked_profile(xxxx, yyyy, xxxx_error, yyyy_error, wed4):

        plt.errorbar(np.array(xxxx),
                     np.array(yyyy), 
                     xerr=xxxx_error,
                     yerr=yyyy_error,
                     linewidth=0, marker='o', markersize=3, alpha=0.95,
                     elinewidth=1, capsize=0, capthick=1,
                     color='black', ecolor='black',
                     label='Stacked image')

        plt.xlabel("Radius, arcmin", fontsize=12)  # "Radius in units of $R_{500}$")
        plt.ylabel("Counts s$^{{-1}}$ arcmin$^{{-2}}$", fontsize=12) # "Brightness in relative units")
        
        plt.xscale("log")
        plt.yscale("log")

        R500inmin = 10

        plt.axvline(R500inmin    , linestyle='--', color='orangered')
        plt.axvline(R500inmin*1.6, linestyle='--', color='dodgerblue')
        plt.axvline(R500inmin*2.7, linestyle='--', color='green')
        plt.axvline(R500inmin*8.1, linestyle='--', color='magenta')
        
#        plt.text(19, 1e-1, '$R_{500c}:R_{200c}:R_{200m}:R_{ta}=$\n$=1:1.6:2.7:8.1$', fontsize=11,
#                 bbox=dict(facecolor='white', alpha=0.99, edgecolor='grey'),
#                 ha='center', va='center')   
                
        plt.legend(loc=3, fontsize=12)
        plt.xticks([0.1, 1, 10, 100], [0.1, 1, 10, 100])
        
        # 4 different plot for 4 wedges:
        if True:
        
            br_all = np.vstack([wed4])
           # print(wed4)
            idx = np.argsort(br_all, axis=0)[::-1, :]
            br_sort = np.take_along_axis(br_all, idx, axis=0)
           # print(br_sort)
  
            br1, br2, br3, br4 = br_sort          
          
            phon = 6e-6
                               
          #  plt.plot(np.array(xxxx), np.array(br1)-phon)
          #  plt.plot(np.array(xxxx), np.array(br2)-phon)
          #  plt.plot(np.array(xxxx), np.array(br3)-phon)
          #  plt.plot(np.array(xxxx), np.array(br4)-phon)
            
            plt.fill_between(np.array(xxxx),
                             np.array(br4)-phon, np.array(br1)-phon,
                             color='magenta', 
                             zorder=10, alpha=0.1)
            plt.fill_between(np.array(xxxx),
                             np.array(br3)-phon, np.array(br2)-phon,
                             color='magenta', 
                             zorder=10, alpha=0.3)           
        
        resc1 = lambda x: x/R500inmin
        resc2 = lambda x: x*R500inmin
        
        ax2 = plt.gca().secondary_xaxis("top", functions=(resc1, resc2))
        ax2.set_xlabel("Radius / R$_{500}$", fontsize=12)
        ax2.set_xticks([0.01, 0.1, 1, 10], [0.01, 0.1, 1, 10])
        
        return None
