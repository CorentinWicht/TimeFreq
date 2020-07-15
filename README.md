# EEG Time-frequency analysis scripts

The MATLAB scripts available in this repository enable to compute EEG time-frequency decomposition as well as statistical analysis of simple designs (up to 1 within-subject and 1 between-subject factors). 



## Getting Started

Start by clicking on `Code` on the top right of the screen and then `Download ZIP` to download the whole repository (alternatively you can also clone it). 
Then, you will need to run the scripts **in the following order**:

```
1.TimeFreq_Design.m
2.TimeFreq_Main.m
3.TimeFreq_Figures.m (optional and still in development)
```

You will find below a step by step explanation of how to run each script.

### 1.TimeFreq_Design.m

This script is the crucial to the proper functioning of the analysis. 

Start by running the script in MATLAB (click on the `Run` button or on `F5`).

X prompts will appear in the following order

#### 1.1. Settings

[![](tools/screenshots/Settings.png)]

The first prompt enables to set up parameters regarding the time-frequency decomposition.

You can decide to compute 1) `Evoked activity`, 2) `Induced activity` and/or 3) `Intertrial coherence`.
For a deeper understanding between these three measures see:
**REFS**

Additionally, you can either import a `TimeFreq_Design.mat` and/or a `TimeFreq_Bands.mat` files (Y)  or create new ones (N).

Finally, you can decide to :
* Average reference your EEG files before the time-frequency decomposition
* Restrict the time-window for statistics inside the time-range of your ERP
* Apply a power spectrum normalization as usually performed by EEGLAB (i.e. 10 x Log10)

## Dependencies


## Author

* **Corentin Wicht** 
*SNSF Doc.CH PhD student*
*[Laboratory for Neurorehabilitation Science](https://www3.unifr.ch/med/spierer/en/), University of Fribourg, Switzerland*
*corentin.wicht@unifr.ch; corentinw.lcns@gmail.com*

## License

This project is licensed under the XXX License - see the [LICENSE.md](LICENSE.md) file for details
 
## Versions

XXX
