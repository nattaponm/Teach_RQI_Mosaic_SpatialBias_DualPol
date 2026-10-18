"""radteach — teaching toolkit for radar QPE, RQI mosaics, dual-pol rainfall and gauge verification.

Course: Teach_RQI_Mosaic_SpatialBias_DualPol (Chiang Mai flood, 3 Oct 2024).
Modules
    config    paths and YAML settings
    colab     Colab setup, data download, saving results to the student's own fork
    io        reading UF radar files, scan inventory
    dem       DEM tile mosaic, hillshade, terrain interpolator
    geometry  beam height, beam radius, beam blockage (CBB)
    qc        quality control, ZDR/PhiDP calibration, KDP, attenuation
    qpe       Z-R, R(KDP), R(A), blended dual-pol
    rqi       radar quality index
    vpr       melting layer and vertical profile of reflectivity
    mosaic    hybrid scans, polar-to-grid, multi-radar mosaic, combine() for MAX/MEAN/NEAREST/RQI_MAX/RQI_WMEAN
    accum     15-min to hourly accumulation, advection correction
    gauges    HII hourly gauges: download, QC, time alignment
    verify    gauge-radar pairing and verification metrics
    bias      mean-field and spatial bias adjustment
    plotting  journal-style figures, Thai fonts, colour maps
    mapping   standard maps (graticule, scale bar, north arrow, inset)
"""
__version__ = "1.0.0"
