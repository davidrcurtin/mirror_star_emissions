import matplotlib.pyplot as plt
import matplotlib as mpl
from module import Nugget, load_thin_data
import numpy as np
import physics

standard_rho_c_MS = [0.01,0.1,1,10,100]
custom_cmap = mpl.colors.ListedColormap(['#0E0880', '#2518F1', '#3352FF', '#54B2FF', '#78FECD', '#9BFE79', '#D5FE28', '#F5C30D', '#ED6912', '#DD1012', '#750205'])
custom_cmap.set_extremes(over= 'k')


def plot_nugget_paper(nugget:Nugget, plotpath = None, figsize=(4.5,6.3)):
    '''Create the plot used in the paper for one nugget: profiles of T, rho, Prad/Pgas. 
    If plotpath is not provided, plot will be shown, otherwise plot will be saved to the given path.'''

    rho_c_MS_frac = nugget.rho_c_MS_frac()
    xi = nugget.xi()
    r = nugget.r(km=True)
    r_photo = nugget.r_photo(km=True)
    T = nugget.T()
    rho = nugget.rho(cgs=True)
    M_photo = nugget.M_photo()
    L_photo = nugget.L_photo()
    L_heating = nugget.L_heat()
    convective_tag = nugget.convective_tag()

    fig, axs = plt.subplots(ncols=1, nrows=3, sharex=True, figsize=figsize, dpi=300, gridspec_kw={'hspace': 0})

    # T plot
    ax = axs[0]
    ax.scatter(r, T, c = -1.0*convective_tag, s= 1, cmap='flag', label = 'Convective')
    ylim = ax.get_ylim()
    ax.scatter(-2,1, s=1, c='k', label = 'Radiative')
    ax.set_ylim(ylim)
    ax.axvline(r_photo, c = 'blue', linestyle = '--', label = 'Photosphere')
    ax.set_ylabel(r'$T\ \rm[K]$')
    ax.set_xlabel(r'$r\ \rm[km]$')

    # Make legend with size of markers all large
    lgnd = ax.legend(loc="upper right")
    for handle in lgnd.legend_handles[:2]:
        handle.set_sizes([20.0])


    # rho plot
    ax = axs[1]
    ax.scatter(r, rho, c = -1.0*convective_tag, s= 1, cmap='flag', label = 'Convective')
    ax.axvline(r_photo, c = 'blue', linestyle = '--', label = 'Photosphere')
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2e'))
    ax.set_ylabel(r'$\rho\ \rm[g/cm^3]$')
    ax.set_xlabel(r'$r\ \rm[km]$')

    # prad/pgas plot
    ax = axs[2]
    ax.scatter(r, physics.P_gas(T, nugget.rho())/nugget.P_rad(), c = -1.0*convective_tag, s= 1, cmap='flag', label = 'Convective')
    ax.axvline(r_photo, c = 'blue', linestyle = '--', label = 'Photosphere')
    ax.set_ylabel(r'$P_{\rm{gas}}/P_{\rm rad}$')
    ax.set_xlabel(r'$r\ \rm[km]$')

    xi_exp = int(np.floor(np.log10(abs(xi))))
    xi_mantissa = xi/10**xi_exp
    M_exp = int(np.floor(np.log10(abs(M_photo))))
    M_mantissa = M_photo/10**M_exp
    L_exp = int(np.floor(np.log10(abs(L_photo/physics.L_sun))))
    L_mantissa = (L_photo/physics.L_sun)/10**L_exp
    xi_title_piece = rf'$\xi = 10^'+r'{'+f'{xi_exp}'+ r'}$, '
    rho_title_piece = r'$\rho_{\rm{core}} = ' + f'{rho_c_MS_frac}' + r'\rho_{\rm{core},\odot}$' +'\n'
    Lphoto_title_piece = r'$\mathcal{L}_{\rm{photo}}='  rf'{L_mantissa:.2f}\cdot 10^'+r'{'+f'{L_exp}'+ r'}\rm{L_\odot}$,'
    Lratio_title_piece = '' #r'$\mathbf{L}=' +rf'{L_photo/L_heating:.4f}$'+'\n' # Not used in paper anymore
    M_title_piece = r'$\rm{M_{nugget}}=' + rf'{M_mantissa:.2f}\cdot 10^'+r'{'+f'{M_exp}'+ r'} \rm{kg}$'
    title = xi_title_piece + rho_title_piece + Lphoto_title_piece + Lratio_title_piece + M_title_piece
    axs[0].set_title(title)
    
    plt.tight_layout()
    if plotpath is not None:
        plt.savefig(plotpath,bbox_inches = "tight")
    else: plt.show()
    plt.close()

def plot_inputspace(nuggets: list[Nugget], plotpath = None, figsize=(5,4), cbar_range=5):
    '''Make a plot of the given nuggets in the input parameter space (rho_c, T_c) for a single (rho_c_MS, xi) pair.
        cbar_range gives the distance from 0 the color bar will be normalized to, values outside this range will be clamped. 
        If plotpath is not provided, plot will be shown, otherwise plot will be saved to the path given.'''
    rho_c_MS_frac = nuggets[0].rho_c_MS_frac()
    xi = nuggets[0].xi()

    rho_cs = np.log10([n.rho_c(cgs=True) for n in nuggets])
    T_cs = np.log10([n.T_c() for n in nuggets])
    L_ratios = np.array([n.log_L_ratio() for n in nuggets])

    plt.figure(dpi=400, figsize=figsize)
    plt.scatter(rho_cs, T_cs, c= L_ratios, s= 20,zorder = 0, cmap='seismic', vmin = -cbar_range, vmax = cbar_range)

    plt.xlabel(r'$\log_{10} \rho_c [g/cm^3]$')
    plt.ylabel(r'$\log_{10} T_c [K]$')
    plt.colorbar(label = r'$\log_{10}(\mathbf{L})$')
    plt.xlim(min(rho_cs), max(rho_cs))
    plt.ylim(min(T_cs),max(T_cs))

    xi_exp = int(np.floor(np.log10(abs(xi))))
    xi_mantissa = xi/10**xi_exp
    xi_title_piece = rf'$\xi = 10^'+r'{'+f'{xi_exp}'+ r'}$, '
    rho_title_piece = r'$\rho_{\rm{core}} = ' + f'{rho_c_MS_frac}' + r'\rho_{\rm{core},\odot}$' +'\n'
    plt.title(xi_title_piece + rho_title_piece)

    plt.tight_layout()
    if plotpath is not None:
        plt.savefig(plotpath,bbox_inches = "tight")

def plot_paramspace(contour_dict: dict, cmaps = [custom_cmap]*6, 
                    cbar_ranges = [(-10, 1), (0, 7+1/3), (0, 7+1/3), (0, 7+1/3), (-18,7+2/3), (0, 7+1/3)], 
                    cbar_ticks = [np.linspace(-18, -3, 6), [0,2,4,6], [0,2,4,6], [0,2,4,6], [-18, -11, -4, 3], [0,2,4,6]],
                    plotpath = None):
    '''
    :contour_dict: dictionary of nugget data as given by mirrorstarmodule.convert_dictionary()
    :cmaps: list of colormaps to be passed to the scatter and colorbar calls in each column ex. ['viridis', 'seismic', colormap_object, ...]
    :cbar_ranges: list of 2-tuples of the endpoints of each colorbar
    :cbar_ticks: list of lists, each specifying the tick locations on a colorbar'''

    def reformat_contour(contour_dict, columns, modifiers):
        ''' Returns a dictionary with keys = rho_c_MS_frac, values = (xi, M, properties) where properties is a list of np arrays with the same length as  xi & M, representing the properties given by columns.
        Columns is a list of indexes for which properties of the contours to include. Note M is in log space along with all the columns by default, and xi is not. \n
        Ex: [1,2,4] to do just rho_c, T_c, and R\n
        currently inds correspond to [xis, rho_cs, T_cs, L_heating, R_photos, T_photos, rho_photos, m_photos, r_rads, L_ratios, R_max, M_max, tau_center, g_photo] \n
        if modifiers is provided, it will apply each element of modifiers to the data in the same index in the columns parameter, doing nothing if the element is None.'''
        nothing = lambda x : x
        modifiers = [nothing if modifier is None else modifier for modifier in modifiers]

        data_dict = {}
        for rho_c_MS_frac, contours in contour_dict.items():
            contour_data = [[] for i in range(len(columns))]
            xis = []
            Ms = []
            for contour in contours:
                xis.extend(contour[0])
                Ms.extend(contour[7])
                data_ind = 0
                for i, column_ind in enumerate(columns):
                    contour_data[data_ind].extend(modifiers[i](contour[column_ind]))
                    data_ind += 1
            contour_data = [np.array(property) for property in contour_data]
            data_dict[rho_c_MS_frac] = (np.array(xis), np.array(Ms), contour_data)
        return data_dict

    modifiers = [None, None, None, lambda x:x-3, lambda x:np.log10((10**x)/physics.L_sun), lambda x:x+2]

    data = reformat_contour(contour_dict, [1,2,5,4,3,13], modifiers=modifiers)

    col_labels = [r'$\log_{10}(\rho_c) [g/cm^3]$', r'$\log_{10}(T_c) [K]$', r'$\log_{10}(T_{photo}) [K]$', r'$\log_{10}(R_{photo}) [km]$', r'$\log_{10}(L_{photo}/L_{sun})$', r'$\log_{10}(g_{photo}) [cm/s^2]$']
    thin_strings=['rho_c', 'T_c', 'T_c', 'R', 'L', 'g']

    thin_data_trim, thin_data_interp = load_thin_data()
    for rho_c_MS in thin_data_trim.keys():
        thin_data_trim[rho_c_MS]['xi'] = np.log10(thin_data_trim[rho_c_MS]['xi'])
        thin_data_trim[rho_c_MS]['M_nugget'] = np.log10(thin_data_trim[rho_c_MS]['M_nugget'])
            
        thin_data_interp[rho_c_MS]['xi'] = np.log10(thin_data_interp[rho_c_MS]['xi'])
        thin_data_interp[rho_c_MS]['M_nugget'] = np.log10(thin_data_interp[rho_c_MS]['M_nugget'])

    rows, cols = len(data.keys()), len(next(iter(data.values()))[2])
    rho_cores = sorted(data.keys())

    fig, axs = plt.subplots(nrows=rows+1, ncols=cols, figsize=(0.2+cols*2.75, 1.45+rows*2), dpi =500,
                            height_ratios= [0.1]+[1 for i in range(rows)], layout="constrained")
    ax_list = np.ravel(axs[1:, :])
    for ax in ax_list[1:]:
        ax.sharex(ax_list[0])
        ax.sharey(ax_list[0])

    # Plot on each subplot
    for i in range(cols):
        norm = mpl.colors.Normalize(*cbar_ranges[i])
        for j, rho_c_MS in enumerate(rho_cores):
            xi, M, plt_data = np.log10(data[rho_c_MS][0]), data[rho_c_MS][1]+3, data[rho_c_MS][2][i]
            ax = axs[j+1][i] #Leaving space for the colorbar

            # Plot Thick Data
            ax.scatter(xi, M, c = plt_data, s = 10, marker = 'o', norm = norm, cmap = cmaps[i])

            # Plot thin data
            if thin_strings[i] is not None and rho_c_MS in standard_rho_c_MS:
                thin_size = 30
                thin_linewidth = 0.2
                if thin_strings[i] == 'Blank':
                    ax.scatter(thin_data_trim[rho_c_MS]['xi'], thin_data_trim[rho_c_MS]['M_nugget'], c = 'gray',
                                alpha = 0.75, marker = '*', linewidths = thin_linewidth, s=thin_size)
                    ax.scatter(thin_data_interp[rho_c_MS]['xi'], thin_data_interp[rho_c_MS]['M_nugget'], c = 'gray',
                                alpha = 0.75, marker = '$\u2606$', linewidths = thin_linewidth, s=thin_size)
                else:
                    thin_plot_data_trim, thin_plot_data_interp = np.log10(thin_data_trim[rho_c_MS][thin_strings[i]]), np.log10(thin_data_interp[rho_c_MS][thin_strings[i]])

                    ax.scatter(thin_data_trim[rho_c_MS]['xi'], thin_data_trim[rho_c_MS]['M_nugget'], c = thin_plot_data_trim,
                                alpha = 0.75, norm = norm, cmap = cmaps[i], marker = '*', linewidths = thin_linewidth, s=thin_size)
                    ax.scatter(thin_data_interp[rho_c_MS]['xi'], thin_data_interp[rho_c_MS]['M_nugget'], c = thin_plot_data_interp,
                                alpha = 0.75, norm = norm, cmap = cmaps[i], marker = '$\u2606$', linewidths = thin_linewidth, s=thin_size)
            
            # Get rid of interior labels
            ax.tick_params(labelbottom=False, labelleft=False)
        
        # Re-add exterior labels
        axs[-1][i].tick_params(labelbottom=True)
        axs[-1][i].set_xticks(np.arange(-26, -14, 2))

        # Add colorbars
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmaps[i]), location = 'top', cax = axs[0][i], extend = 'max', ticks = cbar_ticks[i])
        cbar.set_label(col_labels[i], fontsize='xx-large', labelpad=12)

    # Add rho_c_MS labels
    for j, rho_c_MS in enumerate(rho_cores):
        ax = axs[j+1][0]
        if rho_c_MS in standard_rho_c_MS:
            rho_c_MS_str = f'{rho_c_MS}'
        else:
            rho_c_MS_str = r'$10^{' + f'{int(np.log10(rho_c_MS))}' + r'}$'
        ax.annotate(rho_c_MS_str+r'$\rho_{c,\odot}$', (0.97, 0.88), xycoords='axes fraction', fontsize='x-large', horizontalalignment='right')
        ax.tick_params(labelleft=True)

    fig.supxlabel(r"$\log \xi$", fontsize='xx-large')
    fig.supylabel(r"$\log (M_{\mathrm{nugget}}/g)$", fontsize='xx-large')

    if plotpath is not None:
        plt.savefig(plotpath, pad_inches=0.1, bbox_inches = "tight")