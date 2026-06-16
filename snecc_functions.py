
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path.append('/Users/walu/icecube/bd_ppros/')
from Event import Event
from Photon import Photon
os.chdir('/Users/walu/icecube/energy_constraints/snecc/')
import copy
import h5py
import numpy as np
from scipy.interpolate import interp1d

#to calculate cross section weight
from asteria.interactions import InvBetaPar, InvBetaTab
from snewpy.neutrino import Flavor
import astropy.units as u

from scipy.ndimage import zoom
import matplotlib.colors as colors

def load_events(bulkice_outfile_path):
    master_df = pd.read_csv(bulkice_outfile_path)
    HC = 1241  # eV·nm

    events = []

    for run_id, group in master_df.groupby('run_id'):
        # --- Event-level attributes ---
        enu = group['enu'].iloc[0]  # neutrino energy for this run
        ee = group['ee'].iloc[0]    # positron energy for this run
        x, y, z = group[['vertex_x', 'vertex_y', 'vertex_z']].iloc[0]
        time = 0 #adjusting for initial time of 1 second
        runID = run_id

        # --- Build photon list for this event ---
        photons = []
        for _, row in group.iterrows():
            energy_eV = row['photon_energy']  # assuming photon_energy is already in eV
            wavelength_m = (HC / energy_eV)   # convert nm
            pmt = row['pmt_hit']
            photon_time = row['hit_time']*1e-9-1e-3   #adjusting for initial time of 1 second
            photon = Photon(wavelength=wavelength_m, pmt=pmt, runID=runID, time=photon_time)
            photons.append(photon)

        # --- Create event and attach photons ---
        event = Event(enu=enu, ee=ee, x=x, y=y, z=z, runID=runID, time=time, weight=1.0)
        event.add_photon_list(photons)
        event.calculate_multiplicity()

        events.append(event)
    del master_df
    return events

def inflate_event(events, n):
    """
    Inflate a list of Event objects by replicating each event `n` times.

    Parameters
    ----------
    events : list of Event
        List of Event objects to be inflated.
    n : int
        Number of times each event should be copied.

    Returns
    -------
    inflated_events : list of Event
        List containing n copies of each original event.
    """
    inflated_events = []
    for event in events:
        for i in range(n):
            # Create a deep copy so each copy is independent
            event_copy = copy.deepcopy(event)
            # Optionally modify the runID or time to differentiate if needed
            inflated_events.append(event_copy)
    return inflated_events

def build_distribution(h5_path, mass_hierarchy="nmo", time_bin_width=0.01, energy_bin_width=2.0):
    """
    Build normalized time and energy distributions from a supernova neutrino model HDF5 file.

    Parameters
    ----------
    h5_path : str
        Path to the .h5 file (e.g., '.../Tamborra2014_27.0Msun_dir1_NU_E_BAR_with_osc.h5')
    mass_hierarchy : str, optional
        One of {'no_osc', 'nmo', 'imo'}. Determines which spectra to use.
    time_bin_width : float, optional
        Desired time bin width in seconds.
    energy_bin_width : float, optional
        Desired energy bin width in MeV.

    Returns
    -------
    result : dict
        {
            "time_grid": 1D array of new time centers,
            "energy_grid": 1D array of new energy centers (rebinned for spectra),
            "spectra_rebinned": 2D array [energy, time],
            "flux_rebinned": 1D total flux vs time,
            "time_distribution": normalized 1D PDF over time (sum = 1),
        }
    """

    # ---------------------------
    # 1. Load the model data
    # ---------------------------
    with h5py.File(h5_path, "r") as f:
        group_name = list(f["models"].keys())[0]
        g = f[f"models/{group_name}"]

        time = g["time"][:]              # seconds
        energy = g["energy"][:]          # MeV
        flux = g["flux"][:]              # 1D flux (vs time)
        spectra = g["spectra"][mass_hierarchy][:].T  # shape: (energy, time)

    # ---------------------------
    # 2. Re-bin time
    # ---------------------------
    t_min, t_max = time.min(), time.max()
    time_edges = np.arange(t_min, t_max + time_bin_width, time_bin_width)
    time_centers = 0.5 * (time_edges[:-1] + time_edges[1:])

    # Interpolate spectra onto new time centers
    interp_spec_time = interp1d(time, spectra, kind="linear", axis=1, fill_value="extrapolate")
    spectra_time_rebinned = interp_spec_time(time_centers)  # shape (E, new_T)

    # Also interpolate flux to same time centers
    interp_flux = interp1d(time, flux, kind="linear", fill_value="extrapolate")
    flux_rebinned = interp_flux(time_centers)

    # ---------------------------
    # 3. Re-bin spectra in energy
    # ---------------------------
    E_min, E_max = energy.min(), energy.max()
    E_edges = np.arange(E_min, E_max + energy_bin_width, energy_bin_width)
    E_centers = 0.5 * (E_edges[:-1] + E_edges[1:])

    # Interpolate spectrum to coarser energy bins
    interp_spec_energy = interp1d(energy, spectra_time_rebinned, kind="linear", axis=0, fill_value="extrapolate")
    spectra_rebinned = interp_spec_energy(E_centers)  # shape: (new_E, new_T)

    # ---------------------------
    # 4. Build normalized time distribution
    # ---------------------------
    flux_integrated = np.trapz(spectra_rebinned, E_centers, axis=0)
    time_distribution = flux_integrated / np.sum(flux_integrated)

    # ---------------------------
    # 5. Return results
    # ---------------------------
    return {
        "time_grid": time_centers,
        "energy_grid": E_centers,
        "spectra_rebinned": spectra_rebinned,
        "flux_rebinned": flux_rebinned,
        "time_distribution": time_distribution,
    }
def flux_weight(events, snewpy_distribution, Ngen=600549, distance = 5 * u.kpc, radius = 20 * u.m):
    E_grid = (snewpy_distribution['energy_grid'] * u.MeV).to(u.J).value
    spectra = (snewpy_distribution['spectra_rebinned'] * u.MeV).to(u.J).value
    time_grid = (snewpy_distribution['time_grid'] * u.s).value
    luminosity = ((snewpy_distribution['flux_rebinned'] * 1e51 * (u.erg /u.s)).to(u.J / u.s)).value
    distance = distance.to(u.m).value
    radius = radius.value

    #Calculating avg energy for each time bin
    avg_energy = []
    avg_energy_square = []
    for j in range(spectra.shape[1]):
        num = np.trapezoid(E_grid * spectra[:, j], E_grid)
        denom = np.trapezoid(spectra[:, j], E_grid)
        if np.isnan(num) or np.isnan(denom) or denom == 0:
            avg_energy.append(0.0)
        else:
            avg_energy.append(num / denom)
        avg_energy_square.append(np.trapezoid((E_grid**2) * spectra[:, j], E_grid) / denom if denom != 0 else 0.0)

    avg_energy = np.array(avg_energy)
    avg_energy_square = np.array(avg_energy_square)

    # Sample time from flux_weight vs time distribution
    interp_flux = interp1d(time_grid, luminosity, kind="linear", fill_value="extrapolate")
    flux_vals = interp_flux(time_grid)

    # Make sure no negative or nan values
    flux_vals = np.nan_to_num(flux_vals, nan=0.0)
    flux_vals = np.clip(flux_vals, a_min=0, a_max=None)

    #calculating model average energy
    model_avg_energy = np.trapezoid(avg_energy*flux_vals, time_grid)/np.trapezoid(flux_vals, time_grid)

    #calculating alpha
    avg_alpha = (2*avg_energy**2 - avg_energy_square) / (avg_energy_square - avg_energy**2) 
    avg_alpha = np.nan_to_num(avg_alpha, nan=0.0, posinf=0.0, neginf=0.0)  
    model_alpha = np.trapezoid(avg_alpha*flux_vals, time_grid)/np.trapezoid(flux_vals, time_grid)

    # Divide safely: replace NaN, inf, or div-by-zero with 0
    with np.errstate(divide='ignore', invalid='ignore'):
        energy_averaged_flux = np.divide(flux_vals, avg_energy)
        energy_averaged_flux[~np.isfinite(energy_averaged_flux)] = 0.0  # set nan, inf to 0

    #integrating total flux over time grid
    flux_weight = (1/Ngen) * (radius**2 / (distance ** 2)) * np.trapezoid(energy_averaged_flux, time_grid)

    for event in events:
        event.weight = event.weight * flux_weight
    return events, model_avg_energy*6.242e12, model_alpha #converting to MeV

def eff_vol_weight(events, abs_file_path, absorption_sim, N_modules, max_multiplicity, multiplicity_params):
    """
    Calculate W_eff(m) for multiplicities up to max_multiplicity.
    """
    with open(abs_file_path, 'r') as f:
        absorption_lengths = [float(line.strip()) for line in f if line.strip()]

    weights = {}

    for m in range(1, max_multiplicity + 1):
        m_key = m if m in multiplicity_params else 6
        b, c = multiplicity_params[m_key]

        # Mean effective volume across all absorption lengths
        V_eff_all = [b * (1.0 / L) + c for L in absorption_lengths]
        V_eff_mean = np.mean(V_eff_all)

        # Effective volume at provided absorption length
        V_eff_sim = b * (1.0 / absorption_sim) + c

        # Weight for this multiplicity
        W_eff = N_modules * (V_eff_mean / V_eff_sim)
        weights[m] = W_eff

    #  Apply multiplicative weighting
    for e in events:
        if hasattr(e, 'multiplicity') and e.multiplicity is not None:
            m = e.multiplicity
            W_eff = weights[m] if m in weights else weights[max(weights)]

            
            e.weight = e.weight * W_eff

    return events

def cross_section_weight(events):
    ibd = InvBetaPar()
    l = (40 * u.m).to(u.cm).value
    ice_density = 0.92 #g/cm3
    Na = 6.022e23 #1/mol
    n_interaction = 2
    M= 18 #g/mol
    n_target = ice_density * Na * (n_interaction / M)
    #n_target =  6.023e23 * (2 / 18)
    for event in events:
        xs = ibd.cross_section(Flavor.NU_E_BAR, event.enu * u.MeV)
        event.weight = event.weight * xs.value * n_target * l
    return events

def build_mean_energy_alpha_grid(h5_filename, mass_hierarchy='no_osc'):
    """
    Load mean-energy/alpha spectra from spectra_grid.h5,
    normalize each spectrum, compute mean energy and alpha,
    and return a dict containing ONLY the spectrum 
    corresponding to the requested mass hierarchy.

'
    """

    # Map user-friendly keyword to dataset name
    hierarchy_map = {
        "no_osc": "ispec",
        "nmo": "nmo",
        "imo": "imo",
    }

    if mass_hierarchy not in hierarchy_map:
        raise ValueError(f"Invalid mass_hierarchy '{mass_hierarchy}'. "
                         f"Choose from {list(hierarchy_map.keys())}.")

    tag_to_read = hierarchy_map[mass_hierarchy]

    with h5py.File(h5_filename, "r") as hf:
        energy_grid = hf["E_axis"][:]
        mean_energy_groups = [k for k in hf.keys() if k != "E_axis"]

        # Storage containers
        mean_energies = []
        alpha_params = []
        spectra_dict = {}   # dict-of-dicts: [mean_E][alpha] = spectrum

        # Loop through mean energy bins
        for me_key in mean_energy_groups:
            me_val = float(me_key)
            spectra_dict[me_val] = {}

            alpha_groups = list(hf[me_key].keys())

            # Loop through alpha bins
            for alpha_key in alpha_groups:
                alpha_val = float(alpha_key)
                grp = hf[me_key][alpha_key]

                # Load only the selected spectrum (ispec, nmo, or imo)
                spectrum_raw = grp[tag_to_read][:]

                # Normalize spectrum
                norm = np.trapz(spectrum_raw, energy_grid)
                if norm == 0:
                    spectrum_norm = np.zeros_like(spectrum_raw)
                    mean_e = 0
                    alpha = 0
                else:
                    spectrum_norm = spectrum_raw / norm

                    # mean energy
                    mean_e = np.trapz(energy_grid * spectrum_norm, energy_grid)

                    # RMS + alpha
                    mean_e_sq = np.trapz((energy_grid**2) * spectrum_norm, energy_grid)
                    variance = mean_e_sq - mean_e**2

                    if variance <= 0:
                        alpha = 0
                    else:
                        alpha = (2*mean_e**2 - mean_e_sq) / (mean_e_sq - mean_e**2)

                # Store result
                spectra_dict[me_val][alpha_val] = spectrum_norm

                # Bookkeeping
                mean_energies.append(me_val)
                alpha_params.append(alpha_val)

    return {
        "mean_energies": np.array(mean_energies),
        "alpha_params": np.array(alpha_params),
        "spectra": spectra_dict,     # dict[mean_E][alpha] = spectrum(normed)
        "energy_grid": energy_grid,
        "time_grid": None,           
        "mass_hierarchy": mass_hierarchy,
        "tag_used": tag_to_read      # useful for debugging
    }


def spectrum_weight(events, mean_enu, alpha, alpha_meanE_grid):
    """
    Given a value of mean_enu and alpha,
    it calculates the associated spectrum 
    weight assuming the SN spectrum does 
    not change with time
    """
    energy_grid = alpha_meanE_grid['energy_grid']
    normalized_spectrum = alpha_meanE_grid['spectra'][mean_enu][alpha]
    #weighting the events with the spectrum at the sampled time
    weighted_events = []
    for i, event in enumerate(events):
        #finding closest time to the sampled time
        interp_spec = lambda E: np.interp(E, energy_grid, normalized_spectrum, left=0.0, right=0.0)
        weight = interp_spec(event.enu)
        event.weight = event.weight * weight
        weighted_events.append(event)

    return weighted_events

def spectrum_weight_array(base_enu, base_weight0, energy_grid, normalized_spectrum):
    """
    Vectorized spectrum weighting.
    Returns a new weight array: base_weight0 * w(enu)
    """
    # Interpolate weights all at once
    #print(f"Warning: using vectorized spectrum weight")
    spec_weights = np.interp(base_enu, energy_grid, normalized_spectrum,
                             left=0.0, right=0.0)
    
    if np.any(np.isnan(spec_weights)):
        print(f"WARNING: Nan weight found for base_enu: {base_enu[np.isnan(spec_weights)]} with base_weight: {base_weight0[np.isnan(spec_weights)]}")
        spec_weights = np.nan_to_num(spec_weights, nan=0.0)
    return base_weight0 * spec_weights


def QE_weights(events, qe_file_path="data/qe_data.txt", min_multiplicity = 1, max_multiplicity = 24):
    wavelength, qe = np.loadtxt(qe_file_path, dtype=float, unpack=True)
    wavelength_pdf = qe  #QE should not be normalized, it does not have to be
    interp_qe = interp1d(wavelength, wavelength_pdf, kind="linear", fill_value=0.0, bounds_error=False)
    for event in events:
        photons = event.get_photons()
        p_list = np.zeros(len(photons))
        for i, photon in enumerate(photons):
            p_list[i] = interp_qe(photon.wavelength)
        total_p = np.mean(p_list)
        event.weight = event.weight * total_p  

    return events 

def plot_ratio_grid(meanE, alpha, multiplicities, weights, m_weight, gt_weight, multiplicity):
    multiplicities = np.array(multiplicities)
    weights = np.array(weights)
    ratio = m_weight/gt_weight
    # Histogram with weights
    bins = np.arange(min(multiplicities), max(multiplicities) + 2) - 0.5

    xticks = np.arange(0, 20, 1)             # or whatever ticks you want
    yticks = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6, 1e7]   # for log-scale
    plt.figure(figsize=(7, 5))
    plt.hist(
        multiplicities,
        bins=bins,
        weights=weights,
        edgecolor='black',
        alpha=0.75,
        color='royalblue',
        label=f"{meanE:.2f} ratio: {ratio:.4f}\n gt: {gt_weight:.2e}\n m: {m_weight:.2e}\n multiplicity: {multiplicity}"    
    )

    plt.legend(loc='upper right')
    plt.yscale('log')

    # set limits BEFORE saving
    plt.xlim([0, 15.5])
    #plt.ylim([1e-10, 1e7])

    # set ticks BEFORE saving
    #plt.xticks(xticks)
    #plt.yticks(yticks)
    plt.savefig(f'/Users/walu/icecube/energy_constraints/temp_plots/multiplicity_meanE_{meanE:.2f}_alpha_{alpha:.2f}.png', dpi = 300)
    plt.close()

def calculate_coincidence_ratio(
    events,
    snewpy_distribution,
    multiplicity_min=2,
    multiplicity_max=12,
    filename='spectra_grid.h5',
    ref_filename = 'ref_spectra_grid.h5',
    mass_hierarchy='no_osc',
    plot_ratio=False,
    save_multiplicity=2,
    Ngen=600549, 
    distance = 5, 
    radius = 20,
    N_modules=10000
):
    
    distance = distance * u.kpc
    radius = radius * u.m
    
    #Returns coincidence ratio grids for a range of multiplicities.

    #Output:
     #   ratio_grids[m] -> 2D array R[i,j] for multiplicity m

    # Load the spectra grid
    alpha_meanE_grid = build_mean_energy_alpha_grid(filename)
    mean_energies = np.unique(alpha_meanE_grid['mean_energies'])
    alphas        = np.unique(alpha_meanE_grid['alpha_params'])

    #loading reference spectra grid
    #ref_alpha_meanE_grid = build_mean_energy_alpha_grid(ref_filename)
    #ref_mean_energies = np.unique(ref_alpha_meanE_grid['mean_energies'])
    #ref_alphas        = np.unique(ref_alpha_meanE_grid['alpha_params'])

    # ------------------------------------------------
    # Apply all GLOBAL weights ONCE
    # ------------------------------------------------
    base_events = inflate_event(events, 1)
    base_events = cross_section_weight(base_events)
    base_events, model_avg_energy, model_alpha = flux_weight(
        base_events,
        snewpy_distribution=snewpy_distribution,
        Ngen=Ngen, 
        distance = distance, 
        radius = radius
    )
    #base_events = QE_weights(base_events)
    print(f"WARNING: QE weights are currently not applied.\n")

    multiplicity_params = {
        1: (15.1843, 0.000),
        2: (0.1352, 47.005),
        3: (0.0306, 16.170),
        4: (0.0095, 8.680),
        5: (0.0039, 5.105),
        6: (0.0026, 2.909),
    }

    base_events = eff_vol_weight(
        events=base_events,
        abs_file_path='/Users/walu/icecube/energy_constraints/snecc/data/absorption.txt',
        absorption_sim=227.2,
        N_modules=N_modules,
        max_multiplicity=24,
        multiplicity_params=multiplicity_params
    )
    print(f"WARNING: absorption file path is hardcoded.\n")
    for event in base_events:
        if np.isnan(event.weight):
            event.weight = 0.0  
    """
    #Weighting the reference events
    ref_base_events = inflate_event(events, 1)
    ref_base_events = cross_section_weight(ref_base_events)
    ref_base_events = flux_weight(
        ref_base_events,
        snewpy_distribution=snewpy_distribution,
        Ngen=Ngen
    )
    #ref_base_events = QE_weights(ref_base_events)
    ref_base_events = eff_vol_weight(
        events=ref_base_events,
        abs_file_path='../doumeki_analysis/absorption.txt',
        absorption_sim=227.2,
        N_modules=N_modules,
        max_multiplicity=24,
        multiplicity_params=multiplicity_params
    )
    """
    # ------------------------------------------------
    # Prepare output containers
    # ------------------------------------------------
    multiplicities = range(multiplicity_min, multiplicity_max + 1)

    ratio_grids = {
        m: np.zeros((len(mean_energies), len(alphas)))
        for m in multiplicities
    }
    m_weights = {
        m: np.zeros((len(mean_energies), len(alphas)))
        for m in multiplicities
    }
    gt_weights = {
        m: np.zeros((len(mean_energies), len(alphas)))
        for m in multiplicities
    }
    """
    #doing the same for reference weights

    ref_ratio_grids = {
        m: np.zeros((len(ref_mean_energies), len(ref_alphas)))
        for m in multiplicities
    }
    ref_m_weights = {
        m: np.zeros((len(ref_mean_energies), len(ref_alphas)))
        for m in multiplicities

    }
    ref_gt_weights = {
        m: np.zeros((len(ref_mean_energies), len(ref_alphas)))
        for m in multiplicities
    }

    #looping through ratio grids for reference spectra
    ref_base_events = spectrum_weight(ref_base_events, ref_mean_energies[0], ref_alphas[0], ref_alpha_meanE_grid)
    for m in multiplicities:
        ref_m_weight = sum(
                    e.weight for e in ref_base_events
                    if e.multiplicity == m
                )
        ref_gt_weight = sum(
                    e.weight for e in ref_base_events
                    if (e.multiplicity >= multiplicity_min) & (e.multiplicity <= multiplicity_max)
                )
        ref_m_weights[m][0,0] = ref_m_weight
        ref_gt_weights[m][0,0] = ref_gt_weight
        ref_ratio_grids[m][0,0] = (ref_m_weight / ref_gt_weight) if ref_gt_weight > 0 else 0.0
    
    print(f"Finished Calculation of reference ratio\n")
    """
    # ------------------------------------------------
    # Loop through grid
    # ------------------------------------------------
    # Extract base quantities once
    base_enu = np.array([e.enu for e in base_events])
    base_mult = np.array([e.multiplicity for e in base_events])
    base_weight0 = np.array([e.weight for e in base_events])

    # Precompute multiplicity masks (no need to redo every time)
    gt_mask = (base_mult >= multiplicity_min) & (base_mult <= multiplicity_max)
    masks = {m: (base_mult == m) for m in multiplicities}

    #setting up spectrum and energy grid
    energy_grid = alpha_meanE_grid['energy_grid']
    spectra = alpha_meanE_grid['spectra']

    for i, meanE in enumerate(mean_energies):
        for j, alpha in enumerate(alphas):

            # Get the preloaded normalized spectrum for this (meanE, alpha)
            normalized_spectrum = spectra[meanE][alpha]

            # Compute weights for all events at once
            weights = spectrum_weight_array(
                base_enu, base_weight0, energy_grid, normalized_spectrum
            )

            # Compute total weight for multiplicities >= threshold
            gt_weight = weights[gt_mask].sum()

            for m in multiplicities:
                m_weight = weights[masks[m]].sum()

                ratio = (m_weight / gt_weight) if gt_weight > 0 else 0.0

                ratio_grids[m][i, j] = ratio
                m_weights[m][i, j] = m_weight
                gt_weights[m][i, j] = gt_weight
                if (m == save_multiplicity) & (plot_ratio):
                    plot_ratio_grid(meanE, alpha, base_mult, weights, m_weight, gt_weight, save_multiplicity)


    return {
        "mean_energies": mean_energies,
        "alphas": alphas,
        "ratio_grids": ratio_grids,   # dict: multiplicity → 2D grid
        "m_weight": m_weights,        # dict: multiplicity → 2D grid
        "gt_weight": gt_weights,       # dict: multiplicity → 2D grid
        'model_avg_energy': model_avg_energy,
        'model_alpha': model_alpha
    }

"""
m_weight is number of events with the given multiplicity
gt_weight is the number of events with multiplicity greater 
than the given one
"""

def std_of_ratio(ratio_result, multiplicity, min_multiplicity, ref_meanE, ref_alpha, background_path = '../../OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv', total_time = 13 * 60):
    """
    Calculates the standard deviation of reference coincidence ratio
    currently assumes zero background, will update soon.
    """
    background_data = pd.read_csv(background_path)
    try:
        background_m = background_data[background_data['multiplicity'] == multiplicity]['rate'].values * total_time
        print(f"Background for multiplicity {multiplicity}: {background_m:.2f} events over {total_time/60:.1f} minutes")
    except:
        background_m = 0.0
    #this background will be substituted with realistic bkg later
    background_gt = 0
    for rate in background_data[background_data['multiplicity'] > min_multiplicity]['rate'].values:
        background_gt += rate * total_time
    sigma_bkg_m = np.sqrt(background_m)
    sigma_bkg_gt = np.sqrt(background_gt)

    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    meanE_idx = (np.abs(mean_energies - ref_meanE)).argmin()
    alpha_idx = (np.abs(alphas - ref_alpha)).argmin()

    ref_m = ratio_result['m_weight'][multiplicity][meanE_idx, alpha_idx]
    ref_gt = ratio_result['gt_weight'][multiplicity][meanE_idx, alpha_idx]

    sigma_m = np.sqrt(ref_m + sigma_bkg_m**2)
    sigma_gt = np.sqrt(ref_gt + sigma_bkg_gt**2)
    ref_R = ratio_result['ratio_grids'][multiplicity][meanE_idx, alpha_idx]

    #propagating the errors
    ref_sigma = np.sqrt((sigma_m / ref_gt)**2 + (ref_R**2 * (sigma_gt / ref_gt)**2))

    return ref_sigma, ref_R


def chi_square_grid_calculation(
        ratio_result, 
        multiplicity=2,
        min_multiplicity = 2,  
        background_file_path='../../OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv'):
    """
    Calculates the chi square grid for a given multiplicity
    compared to a reference mean energy and alpha.
    """
    background_data = pd.read_csv(background_file_path)

    # Extract grids and parameters
    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    ratio_grid = ratio_result['ratio_grids'][multiplicity]


    # Reference values
  
    #ref_R = ratio_result['ref_ratio_grid'][multiplicity]
    #ref_meanE = ratio_result['ref_mean_energy']
    #ref_alpha = ratio_result['ref_alphas']
    ref_meanE = ratio_result['model_avg_energy']
    ref_alpha = ratio_result['model_alpha']
    #printing warning that this function needs to be updated
    print("Warning: using the mean of energies and alphas as mean value.")
    ref_sigma, ref_R = std_of_ratio(ratio_result, multiplicity, min_multiplicity=min_multiplicity, ref_meanE=ref_meanE, ref_alpha=ref_alpha)

    # Initialize chi-square grid
    chi_square_grid = np.zeros_like(ratio_grid)

    # Calculate chi-square for each grid point
    for i in range(len(mean_energies)):
        for j in range(len(alphas)):
            R_ij = ratio_grid[i, j]
            total_sigma = np.sqrt(ref_sigma**2)

            # Avoid division by zero
            if ref_sigma > 0:
                chi_square_grid[i, j] = ((R_ij - ref_R) ** 2) / (total_sigma ** 2)
            else:
                chi_square_grid[i, j] = 0.0

    return chi_square_grid  

def save_coincidence_output(result, outname="coincidence_ratio.h5"):
    """
    Save output of calculate_coincidence_ratio to HDF5.
    """
    with h5py.File(outname, "w") as f:

        # Save 1D grids
        f.create_dataset("mean_energies", data=result["mean_energies"])
        f.create_dataset("alphas", data=result["alphas"])
        f.create_dataset("model_avg_energy", data=result["model_avg_energy"])
        f.create_dataset("model_alpha", data=result["model_alpha"])

        # Grouped 2D grids keyed by multiplicity
        ratio_grp = f.create_group("ratio_grids")
        m_w_grp   = f.create_group("m_weight")
        gt_w_grp  = f.create_group("gt_weight")

        for m, grid in result["ratio_grids"].items():
            ratio_grp.create_dataset(str(m), data=grid)

        for m, grid in result["m_weight"].items():
            m_w_grp.create_dataset(str(m), data=grid)

        for m, grid in result["gt_weight"].items():
            gt_w_grp.create_dataset(str(m), data=grid)

"""
Single chi-square confidence intervals
"""
"""
calculating chi-square fluctuation
and percentile for confidence intervals.
For each grid points, we simulate 10000
poisson fluctuation on m and gt. This works
because sum of poisson is also poisson.
"""
def calculate_single_confidence_intervals(ratio_result, multiplicity=3):
    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    ratio_grid = ratio_result['ratio_grids'][multiplicity]
    ms = ratio_result['m_weight'][multiplicity]
    gts = ratio_result['gt_weight'][multiplicity]

    num_simulations = 10000
    chi_squares = np.zeros((len(mean_energies), len(alphas), num_simulations))
    factor = 1.0

    for i in range(len(mean_energies)):
        for j in range(len(alphas)):
            ref_m = ms[i, j] 
            ref_gt = gts[i, j] 
            ref_R = ratio_grid[i, j]

            ref_sigma, ref_R = std_of_ratio(
                ratio_result,
                multiplicity,
                2,
                mean_energies[i],
                alphas[j]
            )
            ref_R *= factor
            for k in range(num_simulations):
                sim_m = np.random.poisson(ref_m)
                sim_gt = np.random.poisson(ref_gt)

                sim_R = factor * sim_m / sim_gt if sim_gt > 0 else 0.0

                if ref_sigma > 0:
                    chi_squares[i, j, k] = ((sim_R - ref_R) ** 2) / (ref_sigma ** 2)
                else:
                    chi_squares[i, j, k] = 0.0

    #flattening chi-square and calculating percentiles
    chi_squares_flat = chi_squares.ravel()

    ci_68 = np.percentile(chi_squares_flat, 68)
    ci_99 = np.percentile(chi_squares_flat, 99)

    return chi_squares_flat, ci_68, ci_99

"""
Function to calculate std for
combined chi-square method.
"""
def combined_std_of_ratio(ratio_result,
                          multiplicity,  
                          ref_meanE, 
                          ref_alpha,
                          min_multiplicity=2, 
                          max_multiplicity=12, 
                          background_path = '../../OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv', 
                          total_time = 13*60):
    background_data = pd.read_csv(background_path)

    #loading reference values
    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    meanE_idx = (np.abs(mean_energies - ref_meanE)).argmin()
    alpha_idx = (np.abs(alphas - ref_alpha)).argmin()
    ref_m = ratio_result['m_weight'][multiplicity][meanE_idx, alpha_idx]
    ref_gt = ratio_result['gt_weight'][multiplicity][meanE_idx, alpha_idx]
    sum_gt = ref_gt
    #print(f"warning: function still uses old ways of defining reference point in (E, alpha)\n")

    #defining sigma parameters
    background_gt = 0
    for rate in background_data[(background_data['multiplicity'] >= min_multiplicity)&(background_data['multiplicity'] <= max_multiplicity)]['rate'].values:
        background_gt += rate * total_time
    sigma_sum = np.sqrt(sum_gt + background_gt)
    try:
        background_ref = background_data[background_data['multiplicity'] == multiplicity]['rate'].values[0] * total_time
    except:
        background_ref = 0.0
    sigma_ref = np.sqrt(ref_m + background_ref)
    ref_R = ratio_result['ratio_grids'][multiplicity][meanE_idx, alpha_idx]
    cov_ref = -2*ref_R*(sigma_ref/sum_gt)**2
    ref_std=np.sqrt((sigma_ref/sum_gt)**2 + ref_R**2*(sigma_sum/sum_gt)**2+cov_ref)

    return ref_R, ref_std

def combined_chi_square_grid(
    ratio_result,
    multiplicity_min,
    multiplicity_max,
    background_file_path = '/Users/walu/icecube/OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv',
):
    """
    Calculate combined chi-square grid over a range of multiplicities.
    """

    # Extract grids and parameters
    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    min_ratio_grid = ratio_result['ratio_grids'][multiplicity_max]
    # Initialize combined chi-square grid
    combined_chi_square = np.zeros_like(min_ratio_grid)
    del min_ratio_grid

    ref_meanE = ratio_result['model_avg_energy']
    #ref_meanE = np.mean(np.array(ratio_result['mean_energies']))
    ref_alpha = ratio_result['model_alpha']
    #ref_meanE = 11.4
    #ref_alpha = 2.4

    # Loop over multiplicities and accumulate chi-square grids
    for m in range(multiplicity_min, multiplicity_max + 1):
        ref_R, ref_std = combined_std_of_ratio(
            ratio_result=ratio_result,
            multiplicity=m,
            ref_meanE=ref_meanE,
            ref_alpha=ref_alpha,
            min_multiplicity=multiplicity_min,
            max_multiplicity=multiplicity_max,
            background_path=background_file_path)
        # Calculate chi-square for each grid point
        ratio_grid = ratio_result['ratio_grids'][m]
        for i in range(len(mean_energies)):
            for j in range(len(alphas)):
                R_ij = ratio_grid[i, j]
                ref_std = np.sqrt(ref_std**2)
                if ref_std > 0:
                    combined_chi_square[i, j] += ((R_ij - ref_R) ** 2) / (ref_std ** 2)
                else:
                    combined_chi_square[i, j] += 0.0
    return combined_chi_square

"""
Calculating combined chi-square interval
"""
def calculate_combined_confidence_intervals(
    ratio_result,
    min_multiplicity,
    max_multiplicity
):
    mean_energies = ratio_result['mean_energies']
    alphas = ratio_result['alphas']
    num_simulations = 10000

    # total number of flattened entries
    n_E = len(mean_energies)
    n_A = len(alphas)

    chi_squares_flat = np.empty(n_E * n_A * num_simulations, dtype=np.float64)
    idx = 0  # flat index pointer

    for i, E in enumerate(mean_energies):
        for j, alpha in enumerate(alphas):

            # initialize summed chi-square for all simulations
            chi_sum = np.zeros(num_simulations, dtype=np.float64)

            for m in range(min_multiplicity, max_multiplicity + 1):
                ratio_grid = ratio_result['ratio_grids'][m]
                ms = ratio_result['m_weight'][m]
                gts = ratio_result['gt_weight'][m]

                ref_m = ms[i, j]
                ref_gt = gts[i, j]
                ref_R = ratio_grid[i, j]

                ref_R, ref_std = combined_std_of_ratio(
                    ratio_result=ratio_result,
                    multiplicity=m,
                    ref_meanE=E,
                    ref_alpha=alpha,
                    min_multiplicity=min_multiplicity,
                    max_multiplicity=max_multiplicity
                )

                if ref_std <= 0:
                    continue

                # vectorized Poisson draws
                sim_m = np.random.poisson(ref_m, size=num_simulations)
                sim_gt = np.random.poisson(ref_gt, size=num_simulations)

                # ratio with zero-protection (same logic as before)
                sim_R = np.zeros(num_simulations)
                mask = sim_gt > 0
                sim_R[mask] = sim_m[mask] / sim_gt[mask]

                chi_sum += ((sim_R - ref_R) ** 2) / (ref_std ** 2)

            # write directly into flattened array
            chi_squares_flat[idx:idx + num_simulations] = chi_sum
            idx += num_simulations

    ci_68 = np.percentile(chi_squares_flat, 68)
    ci_99 = np.percentile(chi_squares_flat, 99)

    return chi_squares_flat, ci_68, ci_99


#Useful functions but not used in weight calculation
#example loading
def load_coincidence_output(filename="coincidence_ratio.h5"):
    """
    Load coincidence ratio output from HDF5 and return
    the same format as calculate_coincidence_ratio().
    """
    with h5py.File(filename, "r") as f:

        mean_energies = f["mean_energies"][:]
        alphas        = f["alphas"][:]

        ratio_grids = {
            int(m): f["ratio_grids"][m][:]
            for m in f["ratio_grids"].keys()
        }

        m_weights = {
            int(m): f["m_weight"][m][:]
            for m in f["m_weight"].keys()
        }

        gt_weights = {
            int(m): f["gt_weight"][m][:]
            for m in f["gt_weight"].keys()
        }
        model_average_energy = f["model_avg_energy"][()]
        model_alpha = f["model_alpha"][()]

    return {
        "mean_energies": mean_energies,
        "alphas": alphas,
        "ratio_grids": ratio_grids,
        "m_weight": m_weights,
        "gt_weight": gt_weights,
        "model_avg_energy": model_average_energy,
        "model_alpha": model_alpha
    }
import h5py


def load_combined_chisquare(h5_filename):
    """
    Load combined chi-square grid and contour values from HDF5 file.

    Parameters
    ----------
    h5_filename : str
        Path to the chi-square HDF5 file.

    Returns
    -------
    dict containing:
        combined_chi_square_grid : ndarray
        ci68 : float
        ci99 : float
    """

    try:
        with h5py.File(h5_filename, "r") as f:

            if "combined_chi_square_grid" not in f:
                raise KeyError("Dataset 'combined_chi_square_grid' not found.")

            if "ci68" not in f or "ci99" not in f:
                raise KeyError("Contour datasets 'ci68' or 'ci99' not found.")

            combined_grid = f["combined_chi_square_grid"][:]
            ci68 = f["ci68"][()]
            ci99 = f["ci99"][()]

        return {
            "combined_chi_square_grid": combined_grid,
            "ci68": ci68,
            "ci99": ci99
        }

    except FileNotFoundError:
        raise FileNotFoundError(f"Chi-square file not found: {h5_filename}")

    except OSError:
        raise OSError(f"Unable to open HDF5 file: {h5_filename}")

    except Exception as e:
        raise RuntimeError(
            f"Failed to load chi-square data from {h5_filename}"
        ) from e

def plot_grid(mean_energies, mean_alphas, value_grid, smoothing_factor=3, ci_68 = 1, ci_99 = 2):
    """
    Plot 2D colormap of Z = chi-square(meanE, alpha)
    with 68% and 99% confidence contours.
    """

    mean_energies = np.array(mean_energies)
    alphas = np.array(mean_alphas)

    # --- Smoothing ---------------------------------------------------------
    value_smooth = zoom(value_grid, smoothing_factor, order=3)
    
    me_smooth = np.linspace(mean_energies.min(),
                            mean_energies.max(),
                            value_smooth.shape[0])
    alpha_smooth = np.linspace(alphas.min(),
                               alphas.max(),
                               value_smooth.shape[1])

    ME, ALPHA = np.meshgrid(me_smooth, alpha_smooth, indexing="ij")

    # --- Plotting ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.pcolormesh(
        me_smooth,
        alpha_smooth,
        value_smooth.T,
        shading="auto",
        cmap="inferno"
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$log (\chi^2)$", fontsize=14)

    # --- Confidence contours ----------------------------------------------
    levels = [np.log10(ci_68), np.log10(ci_99)]  # 68% and 99%
    colors_contour = ["green", "yellow"]

    contours = ax.contour(
        ME,
        ALPHA,
        value_smooth,
        levels=levels,
        colors=colors_contour,
        linewidths=2.5
    )

    # Label contours directly
    fmt = {
        np.log10(ci_68): "68%",
        np.log10(ci_99): "99%"
    }
    ax.clabel(contours, contours.levels, inline=True, fmt=fmt, fontsize=12)

    # --- Labels ------------------------------------------------------------
    ax.set_xlabel("Mean Energy ⟨E⟩ [MeV]", fontsize=14)
    ax.set_ylabel("Alpha Parameter α", fontsize=14)
    ax.set_title(r"$\chi^2$ Grid with Confidence Contours", fontsize=16)

    #adding grid lines
    ax.grid(visible=True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)

    plt.tight_layout()
    plt.show()
    return contours

def sample_time(events, snewpy_distribution):

    E_grid = snewpy_distribution['energy_grid']
    spectra = snewpy_distribution['spectra_rebinned']  # shape: (E, time)
    time_grid = snewpy_distribution['time_grid']
    flux_weight = snewpy_distribution['flux_rebinned']

    #Calculating avg energy for each time bin
    avg_energy = []
    for j in range(spectra.shape[1]):
        num = np.trapz(E_grid * spectra[:, j], E_grid)
        denom = np.trapz(spectra[:, j], E_grid)
        if np.isnan(num) or np.isnan(denom) or denom == 0:
            avg_energy.append(0.0)
        else:
            avg_energy.append(num / denom)
    avg_energy = np.array(avg_energy)

    # Sample time from flux_weight vs time distribution
    interp_flux = interp1d(time_grid, flux_weight, kind="linear", fill_value="extrapolate")
    flux_vals = interp_flux(time_grid)

    # Make sure no negative or nan values
    flux_vals = np.nan_to_num(flux_vals, nan=0.0)
    flux_vals = np.clip(flux_vals, a_min=0, a_max=None)

    # Divide safely: replace NaN, inf, or div-by-zero with 0
    with np.errstate(divide='ignore', invalid='ignore'):
        energy_averaged_flux = np.divide(flux_vals, avg_energy)
        energy_averaged_flux[~np.isfinite(energy_averaged_flux)] = 0.0  # set nan, inf to 0

    
    flux_pdf = energy_averaged_flux / np.sum(energy_averaged_flux)  # normalize to sum = 1

    #loading the sampled time for each event
    sampled_times = np.random.choice(time_grid, size = len(events), p = flux_pdf)
    for i, event in enumerate(events):
        event.time = sampled_times[i]
        for photon in event.get_photons():
             photon.time += sampled_times[i]

    return events


def calculate_multiplicity(
    events,
    snewpy_distribution,
    requested_params,              # NEW: list of (meanE, alpha) from user
    multiplicity_min=2,
    multiplicity_max=12,
    filename='spectra_grid.h5',
    mass_hierarchy='no_osc',
    save_multiplicity=2,
    Ngen=600549, 
    distance=5, 
    radius=20,
    N_modules=10000
):
    import numpy as np
    import astropy.units as u
    import h5py

    distance = distance * u.kpc
    radius = radius * u.m
    
    # ------------------------------------------------
    # Load spectra grid
    # ------------------------------------------------
    alpha_meanE_grid = build_mean_energy_alpha_grid(filename)
    mean_energies = np.unique(alpha_meanE_grid['mean_energies'])
    alphas        = np.unique(alpha_meanE_grid['alpha_params'])

    # For convenience
    energy_grid = alpha_meanE_grid['energy_grid']
    spectra     = alpha_meanE_grid['spectra']

    # ------------------------------------------------
    # Find nearest (meanE, alpha) grid points to user requests
    # ------------------------------------------------
    def find_closest(val, arr):
        arr = np.asarray(arr)
        return arr[np.argmin(np.abs(arr - val))]

    selected_grid_points = []
    for (E_req, alpha_req) in requested_params:
        closest_E     = find_closest(E_req, mean_energies)
        closest_alpha = find_closest(alpha_req, alphas)
        selected_grid_points.append((closest_E, closest_alpha))

    # remove duplicates
    selected_grid_points = list(set(selected_grid_points))

    print("\nSelected grid points for computation:")
    for E, a in selected_grid_points:
        print(f"  meanE={E}, alpha={a}")

    # ------------------------------------------------
    # Build base events with global weights
    # ------------------------------------------------
    base_events = inflate_event(events, 1)
    base_events = cross_section_weight(base_events)
    base_events, _, _ = flux_weight(
        base_events,
        snewpy_distribution=snewpy_distribution,
        Ngen=Ngen,
        distance=distance,
        radius=radius
    )

    # Effective volume weights
    multiplicity_params = {
        1: (15.1843, 0.000),
        2: (0.1352, 47.005),
        3: (0.0306, 16.170),
        4: (0.0095, 8.680),
        5: (0.0039, 5.105),
        6: (0.0026, 2.909),
    }

    base_events = eff_vol_weight(
        events=base_events,
        abs_file_path='/Users/walu/icecube/energy_constraints/snecc/data/absorption.txt',
        absorption_sim=227.2,
        N_modules=N_modules,
        max_multiplicity=24,
        multiplicity_params=multiplicity_params
    )

    # Clean NaN weights
    for event in base_events:
        if np.isnan(event.weight):
            event.weight = 0.0  

    # ------------------------------------------------
    # Convert to numpy arrays for fast looping
    # ------------------------------------------------
    base_enu  = np.array([e.enu        for e in base_events])
    base_mult = np.array([e.multiplicity for e in base_events])
    base_w0   = np.array([e.weight     for e in base_events])

    # ------------------------------------------------
    # Save to HDF5 (create file if missing)
    # ------------------------------------------------
    import os
    outname = "/Users/walu/icecube/energy_constraints/multiplicity/weighted_multiplicity_histograms.h5"

    # Create file if it does not exist
    if not os.path.exists(outname):
        with h5py.File(outname, "w") as hf:
            hf.attrs["created_by"] = "calculate_multiplicity"
        print(f"Created new HDF5 file: {outname}")
    
    for (meanE, alpha) in selected_grid_points:

        # get normalized spectrum
        normalized_spectrum = spectra[meanE][alpha]

        # Compute new weights array
        weights = spectrum_weight_array(
            base_enu, base_w0, energy_grid, normalized_spectrum
        )

        # assign back into event objects
        for i, e in enumerate(base_events):
            e.weight = weights[i]

        multiplicities_arr = np.array([e.multiplicity for e in base_events])
        weights_arr        = np.array([e.weight       for e in base_events])

        # ------------------------------------------------
        # Weighted multiplicity vector
        # ------------------------------------------------
        hist_weights = np.zeros(multiplicity_max + 1)

        for m in range(multiplicity_min, multiplicity_max + 1):
            mask = (multiplicities_arr == m)
            hist_weights[m] = weights_arr[mask].sum()

        # ------------------------------------------------
        # Write dataset into the HDF5 file
        # ------------------------------------------------
        with h5py.File(outname, "a") as hf:

            # Group name like: meanE_12.500/alpha_2.500
            grp_name = f"meanE_{meanE:.3f}/alpha_{alpha:.3f}"

            # If the group already exists, delete it to overwrite cleanly
            if grp_name in hf:
                del hf[grp_name]

            # Create group and dataset
            grp = hf.create_group(grp_name)
            grp.create_dataset("multiplicity_histogram", data=hist_weights)

        print(f"Saved multiplicity histogram: meanE={meanE}, alpha={alpha}")

    print("\nAll requested multiplicity histograms saved.\n")


"""
def plot_grid_with_contour_paths(
    mean_energies,
    mean_alphas,
    value_grid,
    contours_dict,
    smoothing_factor=3,
    contour_colors=None
):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import zoom

    mean_energies = np.array(mean_energies)
    alphas = np.array(mean_alphas)

    # --- Smoothing ---------------------------------------------------------
    value_smooth = zoom(value_grid, smoothing_factor, order=3)

    me_smooth = np.linspace(mean_energies.min(),
                            mean_energies.max(),
                            value_smooth.shape[0])
    alpha_smooth = np.linspace(alphas.min(),
                               alphas.max(),
                               value_smooth.shape[1])

    # --- Plot --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.pcolormesh(
        me_smooth,
        alpha_smooth,
        value_smooth.T,
        shading="auto",
        cmap="inferno"
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$\log_{10}(\chi^2)$", fontsize=14)

    # --- Plot contour paths ------------------------------------------------
    if contour_colors is None:
        contour_colors = ["lime", "yellow", "cyan", "magenta"]

    for i, (p, contour_list) in enumerate(contours_dict.items()):
        color = contour_colors[i % len(contour_colors)]
        label_done = False

        for x, y in contour_list:
            ax.plot(x, y, color=color, linewidth=2.5)

            # label only once per probability
            if not label_done:
                ax.text(
                    x[len(x)//2],
                    y[len(y)//2],
                    f"{int(100*p)}%",
                    color=color,
                    fontsize=12,
                    weight="bold"
                )
                label_done = True

    # --- Labels ------------------------------------------------------------
    ax.set_xlabel("Mean Energy ⟨E⟩ [MeV]", fontsize=14)
    ax.set_ylabel("Alpha Parameter α", fontsize=14)
    ax.set_title(r"$\chi^2$ Grid with HPD Contours", fontsize=16)

    ax.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    plt.show()

plot_grid_with_contour_paths(
    ratio_result['mean_energies'],
    ratio_result['alphas'],
    np.log10(chi_square_grid + 1e-6),
    contours
)


"""
