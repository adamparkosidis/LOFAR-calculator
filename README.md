# LOFAR-calculator
The LOFAR Unified Calculator for Imaging, or LUCI, is a unified web interface that allows users to compute several parameters relevant for planning an interferometric observation with LOFAR. 

The web interface can be accessed at <https://support.astron.nl/luci>.

LUCI is developed and maintained by the SDC Operations group (SDCO) at ASTRON. For comments and/or feature requests, please contact the SDCO group using the [JIRA Helpdesk](https://support.astron.nl/rohelpdesk).

For a given observation and pipeline setup, the web interface allows users to:

* Compute the raw and pre-processed data size,

* Compute theoretical and effective image sensitivity,

* Estimate the pipeline processing time,

* Plot the visibility profiles of the sky sources on any given date,

* Plan multi-beam observations, and

* Export the observation and pipeline setup in a PDF format.


## Setting up an observation
The Observational setup section in LUCI allows users to specify an interferometric observing run with LOFAR. The various fields under this section are:

* **Observation time (in seconds)** – Duration of an observing run (default: 28800 seconds). LUCI throws an error if a negative value is provided as input.

* **No. of core stations (0 - 24)** – (default: 24)

* **No. of remote stations (0 - 14)** – (default: 14)

* **No. of international stations (0 - 14)** – (default: 14)

* **Number of channels per subband** – Frequency resolution of an observing run. Allowed values are 64, 128, and 256 channels per subband (default: 64)

* **Number of subbands** – Number of subbands to record as part of the observing run (default: 488). LUCI thows an error if the number of subbands is greater than 488.

* **Integration time (in seconds)** – Correlator integration time. This input field determines the time resolution of the observing run (default: 1 second). LUCI thows an error if correlator time is smaller than 0.16 seconds.

* **Antenna set** – specifies the mode in which the stations are operated. Allowed values are LBA Outer, LBA Sparse, HBA Dual, and HBA Dual Inner. (default: HBA Dual Inner).

Upon clicking the button *Calculate*, LUCI computes the *theoretical image sensitivity*, *effective image sensitivity*, *mean elevation* of input target and the *raw data size*. The results of this computation is displayed to the user under the Results section. Note that before doing the calculations, LUCI validates the input specified by the user. If an invalid input is detected, LUCI displays an error.

If the pipeline set-up is also selected, then the Results section will get an additional population of the *processed data size* and *pipeline processing time*.

### Specifying targets, calibrators and A-team sources
Besides calculating the data sizes and processing times, LUCI also allows users to check if their target of interest is visible to LOFAR on a specific date. The target of interest can be specified in the Target field under the Target setup section. If the specified source is valid, the Resolve button can be used to obtain the coordinates of that source. Note that LUCI can resolve valid LoTSS pointing names as well (like P17 or P205\+67). If however, the user wishes to *manually* input the coordinates, a *target name* **must** also be specified.

In addition to plotting the target, users can also plot the elevations of standard LOFAR calibrators and A-team sources by selecting them using the Calibrators and the A-team sources dropdown boxes (these will not be used for data size or processing time calculations).

Upon clicking the *Calculate* button, LUCI produces three plots:

* **Target visibility plot** shows the elevation of a target as seen by LOFAR on the specified date as an interactive plot. If the field No. of international stations (0-14) is set to 0 (i.e.) observing with only the Dutch array, the target elevation is calculated with respect to the LOFAR core. On the other hand, if the user specifies the full array, LUCI plots the minimum apparent elevation of the target as seen by all LOFAR stations. In addition to the user-specified target, LUCI also plots the elevations of the Sun, the Moon, and Jupiter by default. The two blue regions indicate the sunrise and sunset times. A thicker line on top of the input target elevation plot is also shown to highlight the observation time around transit for each target.

* **Beam Layout** plots the tile and station beams for the specified observation.

* **Angular distance** between specified target and other bright sources are presented in a plotly table.

## Exporting the set-up

Users can export their “observational setup” to a PDF file using the *Generate PDF* button. Upon clicking the *Generate PDF* button, LUCI exposes the *Download* file link below the two buttons which can be used to download the generated PDF file.

**Note** that if you click on the *Generate PDF* button before using the *Calculate* button, LUCI will throw an error.

## Installing and running locally
To run the calculator on your local machine, install the dependencies with:
```
pip install -r requirements.txt
```

To launch the calculator, run the following command:
```
python calculator.py
```
Or if the system has python 3, run the following command:
```
python3 calculator.py
```
and point your browser at:
```
0.0.0.0:8051/luci/
```
