# Data processing

## TOC
----
* STEP 1 I0
* STEP 2 BETA
* STEP 3 HOST - SHAPE
* STEP 4 FCIDE - TS
* STEP 5 FCIDE - SHAPE
* STEP 6 BAYES - I0/BETA
* STEP 7 YR

----
## CSVs exported

#### STEP 1 I0
* I0/I0GeneratorForFitting.R gives us fixed median I0 in "CSVs/I0/I0_Calculated.csv"

#### STEP 2 BETA
* host
  * Host/HostCurveGenerator.R gives us F in "Fitting/CSVs/Host/{cultivar}/YearlyWorstSevs.csv"
  * from which can find beta for fixed I0

#### fungicide
* F from Frank paper
* from which can find beta for fixed I0


#### STEP 3 HOST - SHAPE
* Host/HostCurveGenerator.R gives us all control info in "Fitting/CSVs/Host/{cultivar}/FullFrame.csv"
  * this used to fit shape params by fit_model.py in Run

#### STEPS 4 & 5 - FCIDE
* use median I0
* beta from step 2

###### Step 4- TS
* Currently just setting at 240 :(

###### Step 5 - SHAPE
* Extract from "R_input/FrankData.csv"
* write to "CSVs/Fungicide/FrankControl.csv"
* this used to fit shape params by fit_model.py in Run


----
*setting up running model*

#### STEP 6 BAYES - I0/BETA

Bayesian/posterior_fns.py & run_bf.py

* takes in "I0/I0_calculated.csv" for prior
* takes in "Locations/{location}/WorstCults_TopProp6.csv" for severities
* uses these to find betas and I0s, sends to: "Locations/{location}/Generated{paramname}_{setup}.csv"

#### STEP 7 YR
* takes in "R_input/sonderborg_yield_vs_disease_severity"
* fits LM, exports data and line params to 2 files in "CSVs/YR/"


----
## CSVs exported

AllHostData/allData.csv
* by Host/HostCurveGenerator.R, which calls Host/functions.R
* used in Locations/exportWorstCultsSingleLoc.R
* used in Locations/choiceOfLoc.R

Fungicide/FrankControl.csv
* by Fungicide/FungicideControlCurveGenerator.R
* used in fit_model.py to find fcide shape params

Host/{cultivar}/YearlyWorstSevs.csv
* by Host/HostCurveGenerator.R, which calls Host/functions.R
* used to find betas for host fitting

Host/{cultivar}/FullFrame.csv
* by Host/HostCurveGenerator.R, which calls Host/functions.R
* used in fit_model.py to find host shape params

Host/whichHost/MoreThan8Years.csv
* by choiceOfHost
* informs which cultivars we choose

I0/I0_Calculated.csv
* by I0/I0GeneratorForFitting.R
* used in fit_model when finding shape params and betas for fitting

Locations/{locations}/WorstCults_TopProp6.csv
* by Locations files
* used by Bayesian python stuff to find:

Locations/{locations}/Generated_{parname}_{setup}.csv
* used when running model

Locations/whichLoc/MoreThan1000Points.csv
* by choiceOfLoc.R
* informs which locations we use

YR/YieldRel.csv
* YR/linearModelPars.csv
  * by YieldRel/yield_relationship.R
  * informs the link when running model between DS and yield
