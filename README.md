# DCN Modelling

Dorsal cluster neurons axon targeting modelling

The code in this repository implements the data analysis, modelling and simulation presented in "A sequence of probabilistic processes determines differential axon targeting" by Andriatsilavo et al., 2024.

Experimental data files: 
- Fig4E_Fibers-time_TubGFP_Live_updated.xlsx: file that contains the experimental data needed to run the scripts `analyse_experimental_data.py`, `estimate_birth_death_parameters.py`, `simulate_medulla_targeting.py` and `simulate_medulla_targeting_alternativeVersion.py`  
- DATA-set_Analyse_Intravital.xlsx: file that contains the experimental data needed to run the scripts `data_frame_rate_analysis_intraVital_dataset1.py`, `data_frame_rate_analysis_intraVital_dataset2.py`, `data_lifespan_analysis_intraVital_dataset1.py` and `data_lifespan_analysis_intraVital_dataset2.py`   
- dataset-20231211.csv: file that contains the experimental data needed to run the scripts `data_frame_rate_analysis_allTimes_exVivo.py` and `data_lifespan_analysis_exVivo.py`



Code files: 
- `analyse_experimental_data.py`: processes and analyses the experimental data, plots the distributions for the numbers of filopodia and the theoretical and experimental probability of axon axon survival given its filopodia numbers, used for testing for filopodia independence in an axon
- `estimate_birth_death_parameters.py`: calculates the $\kappa$ and $\delta$ parameters of the initially assumed simple birth-death model for describing the filopodia extension and retraction
- `simulate_medulla_targeting.py`: implements the model for filopodium dynamics with $f(t)=\frac{t}{1+e^{-(t-4)}}$(see Supplementary Mathematical Methods), and performs simulations as to reproduce the experimental data. Includes the simulations' results processing and analysis.
- `simulate_medulla_targeting_alternativeVersion.py`: alternative verion of `simulate_medulla_targeting.py` with $f(t)=t$ (see Supplementary Mathematical Methods)
- `data_frame_rate_analysis_allTimes_exVivo.py`: performs the analysis for benchmarking the frame rate used for the experimental acquisition, for the experimental data contained in the data file dataset-20231211.csv
- `data_frame_rate_analysis_intraVital_dataset1.py`: performs the analysis for benchmarking the frame rate used for the experimental acquisition, for the experimental data contained in the tab "RAW_20170817_Ctrl" of the data file DATA-set_Analyse_Intravital.xlsx
- `data_frame_rate_analysis_intraVital_dataset2.py`: performs the analysis for benchmarking the frame rate used for the experimental acquisition, for the experimental data contained in the tab "RAW_20181024_Ctrl" of the data file DATA-set_Analyse_Intravital.xlsx
- `data_lifespan_analysis_exVivo.py`: calculates and plots the filopofia lifespan trajectories for the experimental data contained in the data file dataset-20231211.csv
- `data_lifespan_analysis_intraVital_dataset1.py`: calculates and plots the filopofia lifespan trajectories for the experimental data contained in the tab "RAW_20170817_Ctrl" of the data file DATA-set_Analyse_Intravital.xlsx
- `data_lifespan_analysis_intraVital_dataset2.py`: calculates and plots the filopofia lifespan trajectories for the experimental data contained in the tab "RAW_20181024_Ctrl" of the data file DATA-set_Analyse_Intravital.xlsx

Required modules: os, pandas, numpy, scipy, matplotlib, copy, math








