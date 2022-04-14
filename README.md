# shear_calibration

Reliability-based calibration of Eurocode and _fib_ Model Code shear resistance formulas.

<img src="https://drive.google.com/uc?export=view&id=1FSc6EjO-gDQ_jwg3CR89--X43iJCHAFR" alt="beta_plot" width="500"/>
https://drive.google.com/file/d//view?usp=sharing

- Supporting code to the paper: Slobbe A., Rozsas A., Yuguang Y. () A reliability-based calibration of shear resistance formulas for reinforced concrete members without shear reinforcement (under review).
- Mixing Matlab, Python, and R due to time constraints.
- If you have a question related to the code please open an issue.


## Dependencies

Matlab:
 - developed under Matlab 2021b (earlier under 2018b)
 	* Statistics and Machine Learning Toolbox
 	* Optimization Toolbox
 	* Global Optimization Toolbox
 	* Parallel Computing Toolbox (to reduce wall clock time)
 	* (only for testing: Deep Learning Toolbox (`combvec`))
 - add to your path:
 	* [custom_FERUM](https://github.com/TNO/custom_FERUM): the entire repository
 	* [Statistics---Matlab](https://github.com/rozsasarpi/Statistics---Matlab): the content of the `distribution_functions\univariate\` folder
 	* [Plotting-Matlab](https://github.com/rozsasarpi/Plotting-Matlab): the content of the `distribution_functions\univariate\` folder (used to make prettier plots, non-essential)
 	* [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)


----

## Overview

The header names follow that of the folders.

### `pre_processing` 

#### `model_uncertainty`

* MatLab 
* simple statistical analysis to infer model uncertainty based on experiments

#### `random_variable_algebra`

* Python
* generate pdf and cdf of the product of random variables


### `calibration` 

* MatLab 
* reliability-based calibration
  - partial factor(s): `main_calibration.m`
  - representative value in semi-probabilistic format: `main_calibration_Ck.m`
* uses results from `pre_processing`

### `post_processing`

* R (because MatLab is inconvenient for plotting)
* visualization of the results + shiny web server for quick and interactive data exploration
* A live version (not necessary the most up to date) of the visualization [webserver is available from here](https://rozsasarpi.shinyapps.io/visualize_calibration_results/).

Additional information may be found in dedicated README files in the respective folders.


## On using the repository

* Running the code: 
    - Matlab: run from the folder of the particular file.
    - R: run from the folder of the particular file (programatically ensured if RStudi is used).
    - python: to be run with the working directory set as the root directory of the repository: `\shear_calibration\`.
* Developed and tested under
	- Windows 10.
	- Python 3.x (`pacal`)
	- R 3.x
	- Matlab 2018a
* Install Python dependencies using the `requirements.txt` file.


If action (wind and snow) random variable inputs change:

1) compute product distributions using `pacal` (`code\pre_processing\random_variable_algebra\pacal_product.py`); output: `*.txt` files
2) Run `code\pre_processing\random_variable_algebra\prepare_pacal_for_ferum.m`; output: `\code\calibration\tmp\*.txt`
3) To reduce runtime copy of the content of the txt files from 2) into the relevant `*_pdf.m` and `*_cdf.m` files in  `\code\calibration\calibration_utils\reliability_analysis\`


## Notes on the probabilistic models

<details>
  <summary>Click to expand!</summary>
  
  ## General notes
  to be added: pacal, etc.
  
  ## Particular models
  * traffic load:
    - model uncertainty (`theta_T`): 
      * based on table 10.3 of [^steenbergen2012], considering all components but the load effect component because that we model separately (`theta_E`)
      * its mean is set to 1.0 because it is just a scaler, we scale the model/load during inverse design
      * our expert judgement: the representative value is assumed to be equal to the mean
  	- time-dependent compoponent (`T`)
  	  * Gumbel, CV; based on [^steenbergen2012]
  	  * the characteristic value of `theta_T*T` has a 1-1/1000 non-exceedance probability after [^ec_traffic]
  	  * from the previous points the `P_repr` value of `T` can be computed (see `pacal_product.py`):
  	    - `x_repr = F^{-1}_{theta_T*T}(1-1/100)`
  	    - `P_repr_T = F_{T}(x_repr)`

  [^steenbergen2012]: Steenbergen, R. D. J. M., Morales Napoles, O., Vrouwenvelder, A.C.W.M. (2012). Algemene veiligheidsbeschouwing en modellering van wegverkeerbelasting voor brugconstructies. Retrieved from Delft

  [^ec_traffic]: CEN. (2003a). Eurocode 1: Actions on structures - Part 2: Traffic loads on bridges.

  ## Improvement ideas

  * improve documentaton
  * implement everything in one language (probably python would be the best choice)

</details>







