try:
    import bottleneck as bn
    nmedian = bn.nanmedian
    nmean = bn.nanmean
    nstd = bn.nanstd
    nmin = bn.nanmin
    nmax = bn.nanmax
    nvar = bn.nanvar
except NameError:
    nmedian = np.nanmedian
    nmean = np.nanmean
    nstd = np.nanstd
    nmin = np.nanmin
    nmax = np.nanmax
    nvar = np.nanvar
