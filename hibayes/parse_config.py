"""
# parse_config.py

Functions for parsing a .ini formatted configuration file

"""

from ConfigParser import SafeConfigParser
import os

#-------------------------------------------------------------------------------

def parse_config(filename):
    """ Open a config.ini file and parse it into a dictionary of config values

    Args:
        filename (str): Filename of config file, including path.

    Returns:
        confdict (dict): Python dictionary of runtime parameters, as generated from config file.
    """

    config = SafeConfigParser()
    config.read(filename)

    confdict = {}

    # -------------------------------------------------------------------------------

    # outdir='chains_o3x3'
    #nc_fit=7
    confdict["nc_fit"]  = config.getint("fit", "n_poly")
    confdict["outdir"]  = config.get("file", "outdir")
    #confdict["logfile"] = config.get("file", "logfile")
    #confdict["comment"] = config.get("misc", "comment")
    confdict["nu_1"] = config.getfloat("misc", "nu_ref")

    #-------------------------------------------------------------------------------

    # Data/simulation parameters
    if config.getboolean("simulation", "simulate_sky"):
        confdict["ledaFreqs"] = None
        confdict["ledaSpec"] = None
        confdict["spectrum_errors"] = None
        confdict["coeffs"] = [float(i) for i in config.get("simulation","coeffs").split(",")]
        confdict["ncoeffs"] = len(confdict["coeffs"])
    else:
        confdict["ledaFreqs"] = config.get("file", "frequency")
        confdict["ledaSpec"]  = config.get("file", "spectrum")
        confdict["spectrum_errors"] = config.get("file", "spectrum_errors")

    confdict["seed_SIM"]      = config.getint("simulation", "seed")
    confdict["A_HI_TRUE"]     = config.getfloat("simulation", "A_HI_TRUE")
    confdict["NU_HI_TRUE"]    = config.getfloat("simulation", "NU_HI_TRUE")
    confdict["SIGMA_HI_TRUE"] = config.getfloat("simulation", "SIGMA_HI_TRUE")
    confdict["plot_truth"]      = config.getboolean("simulation", "plot_truth")

    # Observing parameters

    confdict["tObs"] = config.getfloat("observation", "t_obs")    # Obs. time / sec
    confdict["BW"]   = config.getfloat("observation", "chan_bw")  # Delta Freq. / Hz
    confdict["FREQ_MIN"] = config.getfloat("observation", "freq_min")  # MHz
    confdict["FREQ_MAX"] = config.getfloat("observation", "freq_max")  # MHz

    #-------------------------------------------------------------------------------

    # Priors

    confdict["A_HI_PRIOR"]    = config.get("priors", "A_HI_PRIOR")
    confdict["A_HI_MIN"] = config.getfloat("priors", "A_HI_MIN")
    confdict["A_HI_MAX"] = config.getfloat("priors", "A_HI_MAX")
    confdict["NU_MIN"] = config.getfloat("priors", "NU_MIN")
    confdict["NU_MAX"] = config.getfloat("priors", "NU_MAX")
    confdict["SIGMA_HI_MIN"] = config.getfloat("priors", "SIGMA_HI_MIN")
    confdict["SIGMA_HI_MAX"] = config.getfloat("priors", "SIGMA_HI_MAX")
    confdict["BP_PRIOR_RANGE"] = config.getfloat("priors", "BP_PRIOR_RANGE")

    #-------------------------------------------------------------------------------
    nc_fit = confdict["nc_fit"]
    ncoeffs = confdict["ncoeffs"]
    c = confdict["coeffs"]

    # Set up joint parameters
    confdict["parameters"] = ['A_HI', 'NU_HI', 'SIGMA_HI'] + ['p%i' % ic for ic in range(nc_fit)]
    confdict["plotRanges"] = {'A_HI': [confdict["A_HI_MIN"], confdict["A_HI_MAX"]],
                              'NU_HI': [confdict["FREQ_MIN"], confdict["FREQ_MAX"]],
                              'SIGMA_HI': [0.0, confdict["SIGMA_HI_MAX"]]}
    for ic in range(nc_fit):
        confdict["plotRanges"]['p%i' % ic] = [-confdict["BP_PRIOR_RANGE"], confdict["BP_PRIOR_RANGE"]]

    if confdict["plot_truth"]:
        confdict["plotTruth"] = {'A_HI': confdict["A_HI_TRUE"],
                             'NU_HI': confdict["NU_HI_TRUE"],
                             'SIGMA_HI': confdict["SIGMA_HI_TRUE"]}

        for ic in range(ncoeffs):
            confdict["plotTruth"]['p%i' % ic] = c[ic]
        if nc_fit > ncoeffs:
            for ic in range(ncoeffs, nc_fit):
                confdict["plotTruth"]['p%i' % ic] = -1000.0
    else:
        confdict["plotTruth"]=None

    # Plotting
    confdict["labelDict"] = {'A_HI': r'$A_{HI}/\mathrm{mK}$',
                             'NU_HI': r'$\nu_{HI}/\mathrm{MHz}$', \
                             'SIGMA_HI': r'$\sigma_{HI}/\mathrm{MHz}$'}
    for ic in range(nc_fit):
        confdict["labelDict"]['p%i' % ic] = r'$\log_{10}\left(p_{%s}/\mathrm{K}\right)$' % ic

    #-------------------------------------------------------------------------------

    # Set some MultiNEST parameters here
    confdict["n_live_points"] = config.getint("multinest", "n_live_points")
    confdict["max_modes"] = config.getint("multinest", "max_modes")
    confdict["seed"] = config.getint("multinest", "seed")
    confdict["max_iter"] = config.getint("multinest", "max_iter")
    confdict["evidence_tolerance"] = config.getfloat("multinest", "evidence_tolerance")
    confdict["mode_tolerance"] = config.getfloat("multinest", "mode_tolerance")
    confdict["do_ins"] = config.getboolean("multinest", "do_ins")
    confdict["multimodal"] = config.getboolean("multinest", "multimodal")
    confdict["outstem"] = config.get("multinest", "outstem")

    #-------------------------------------------------------------------------------
    outdir = confdict["outdir"]
    confdict["outputfiles_basename"] = os.path.join(outdir, config.get("multinest", "outstem"))

    return confdict
