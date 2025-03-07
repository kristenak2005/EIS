def get_intensity(line,file,fitted_lines,output_location):
    fit_res = fit_data(file,fitted_lines,line,'int',output_location) # Get fitdata
    m = fit_res.get_map(fitted_lines[f'{line}'][1],measurement='intensity') # From fitdata get map
    return m, fit_res

    save_dir = os.path.join(output_location, 'intensity_files')
    os.makedirs(save_dir, exist_ok=True)
  
    timestamp = m.date.strftime("%Y%m%d_%H%M%S")
    output_filename = os.path.join(save_dir, f'eis_{timestamp}_intensity_{line}.fits')
    m.save(output_filename, overwrite=True)

    return m, fit_res




def plot_eis_fits(line, int, vel, wid, aia_map, output_location, fitted_lines):
    figs = (12,5)
    wid_rat = [1,1,1,1]
    asp = 1/4
    fig = plt.figure(constrained_layout=True, figsize=figs)
    gs = gridspec.GridSpec(1,4,width_ratios=wid_rat)
    plt.rcParams['font.size'] = '10'
    date = int.date.strftime("%Y%m%d_%H%M%S")
    alpha = 0.1
    plot_name = fitted_lines[f'{line}'][2]+'; '+int.date.strftime("%Y/%m/%dT%H:%M:%S")

# Intensity
    ax1 = fig.add_subplot(gs[0,0], projection = int, label = 'a)')

    lwr_bnd = np.percentile(int.data, alpha)
    upr_bnd = np.percentile(int.data, 100-alpha)
    norm = colors.Normalize(vmin=lwr_bnd, vmax=upr_bnd)
    int.plot_settings['norm'] = norm

    int.plot(axes=ax1, title = 'a) Peak intensity', aspect=asp)
    x = ax1.coords[0]
    x.set_axislabel(' ')
    x.set_ticklabel(exclude_overlapping=True)

    plt.colorbar(ax=ax1,location='right', label='')

    # Save the intensity map as a FITS file (optional, if you need it in FITS format too)
    intensity_filename = output_location + '/intensity_files/eis_' + wvl + '_intensity.fits'
    i_map.export_fits(intensity_filename)

