*===============================================================================================
*
* Export into CSV files:
*
* 1. Monthly stock returns
* 2. FF5 + UMD factors/fffactors
*
* Compute residuals from:
*
* 1. CAPM
* 2. FF3
* 3. FF5
* 4. FF5+UMD
*
* Save into a permno-yyyymm file sums of the alphas, fitted, and residuals 
* values for t-12 to t-2
*===============================================================================================

use "CRSPmonthly2019dec_edited.dta", clear			

* drop unnecessary items
keep permno date shrcd exchcd dlstcd dlret ret 

count

gen n=1
collapse (sum) n, by(permno date shrcd exchcd dlstcd dlret ret)
drop n

count

duplicates report permno date

gen long yyyymm = 100*year(date)+month(date)
drop date

* Adjust returns for delisting returns
replace ret    = . if ret<-1

* use -30% for missing delisting returns; except use -55% for Nasdaq (Shumway and Warther)
replace dlret  = -0.3  if missing(dlret) & (dlstcd==500 | (dlstcd>=520 & dlstcd<=584)) & exchcd~=3
replace dlret  = -0.55 if missing(dlret) & (dlstcd==500 | (dlstcd>=520 & dlstcd<=584)) & exchcd==3

replace ret    = (1 + ret) * (1 + dlret) - 1 if ~missing(dlret) & ~missing(ret) 
replace ret    = dlret if ~missing(dlret) & missing(ret) 

drop dlret dlstcd

keep if exchcd>=1 & exchcd<=3 & (shrcd==10 | shrcd==11)
drop exchcd shrcd

drop if missing(ret)

keep if yyyymm>=196307

merge m:1 yyyymm using "fffactors.dta", keep(3) keepusing(rf) nogenerate

gen double retrf = ret - rf
drop ret rf

compress

order permno yyyymm retrf

sort permno yyyymm

export delimited using "returns_monthly.csv", delimiter(",") replace


** FF5 + UDM

use "fffactors.dta", clear

keep yyyymm mktrf smb hml rmw cma umd
keep if yyyymm>=196307 & yyyymm<=201912

sort yyyymm

export delimited using "ff6returns_monthly.csv", delimiter(",") replace



** Read estimates back from the CSV file

import delimited using "estimatedresiduals.csv", delimiter(",") clear

rename v1 yyyymm
rename v2 permno
rename v3 capm_res
rename v4 capm_fit
rename v5 ff3_res
rename v6 ff3_fit
rename v7 ff5_res
rename v8 ff5_fit
rename v9 avgret

rename v10 capm_mktbeta

rename v11 ff3_mktbeta
rename v12 ff3_smbbeta
rename v13 ff3_hmlbeta

rename v14 ff5_mktbeta
rename v15 ff5_smbbeta
rename v16 ff5_hmlbeta
rename v17 ff5_rmwbeta
rename v18 ff5_cmabeta

compress

sort permno yyyymm

save "estimatedresiduals.dta", replace



use "merged_crsp_compustat.dta", clear
keep permno yyyymm retnm exchcd me r12_2 prc

keep if yyyymm>=196306

merge 1:1 permno yyyymm using "estimatedresiduals.dta", nogenerate

keep if ~missing(avgret)

save "mergedresidualdata.dta", replace


* Construct factors	
*forvalues model=0/12 {
forvalues model=0/3 {
	
	disp("Model: `model'")
				
	use "mergedresidualdata.dta", clear

	if `model'==0 {
		rename avgret sortvar
		local dvar = "avgret"
		}
	else if `model'==1 {
		rename capm_res sortvar
		}
	else if `model'==2 {
		rename ff3_res sortvar
		}
	else if `model'==3 {
		rename ff5_res sortvar
		}
	else if `model'==4 {
		rename capm_mktbeta sortvar
		}
	else if `model'==5 {
		rename ff3_mktbeta sortvar
		}
	else if `model'==6 {
		rename ff3_smbbeta sortvar
		}
	else if `model'==7 {
		rename ff3_hmlbeta sortvar
		}
	else if `model'==8 {
		rename ff5_mktbeta sortvar
		}
	else if `model'==9 {
		rename ff5_smbbeta sortvar
		}
	else if `model'==10 {
		rename ff5_hmlbeta sortvar
		}
	else if `model'==11 {
		rename ff5_rmwbeta sortvar
		}
	else if `model'==12 {
		rename ff5_cmabeta sortvar
		}
	else if `model'==13 {
		rename ff5_fit sortvar
		}
	
	keep permno yyyymm me exchcd retnm sortvar prc ff5_???beta
	
	qui keep if ~missing(retnm,me,sortvar)
	
	* compute NYSE breakpoints for the main sorting variable
	qui gen double nyse_sortvar = sortvar if exchcd==1 

	qui bysort yyyymm: egen sv30 = pctile(nyse_sortvar), p(30)
	qui bysort yyyymm: egen sv70 = pctile(nyse_sortvar), p(70)
	
	qui gen byte Qsortvar = 1 if                  sortvar <= sv30 & ~missing(sortvar)
	qui replace  Qsortvar = 2 if sortvar > sv30 & sortvar <= sv70 & ~missing(sortvar)
	qui replace  Qsortvar = 3 if sortvar > sv70                   & ~missing(sortvar)

	* compute small/big breakpoints & 10th percentile
	qui gen double nyse_me = me if exchcd==1 & ~missing(sortvar)
	qui bysort yyyymm: egen me50 = pctile(nyse_me), p(50)
	qui bysort yyyymm: egen me10 = pctile(nyse_me), p(10)
	
	qui gen byte mQsortvar = 1 if me <= me50 & ~missing(me) & ~missing(sortvar)
	qui replace  mQsortvar = 2 if me >  me50 & ~missing(me) & ~missing(sortvar)

	keep permno yyyymm Q* mQ* me retnm me10 prc sortvar ff5_???beta

	qui gen double me_x_retnm = me * retnm

	foreach var of varlist ff5_???beta {
		qui gen me_x_`var' = me * `var'
		}
		
	* UMD-style factor
		
	collapse (sum) me me_x_* (count) Nfirms=me, by(yyyymm mQsortvar Qsortvar)

	qui gen double vwret      = me_x_retnm / me 	
	foreach beta in "mkt" "smb" "hml" "rmw" "cma" {
		qui gen double `beta'beta = me_x_ff5_`beta'beta / me
		}
	
	qui drop if Nfirms==0 

	qui gen     FFport = 1 if mQsortvar==1 & Qsortvar==1	
	qui replace FFport = 2 if mQsortvar==1 & Qsortvar==2	
	qui replace FFport = 3 if mQsortvar==1 & Qsortvar==3	
	qui replace FFport = 4 if mQsortvar==2 & Qsortvar==1	
	qui replace FFport = 5 if mQsortvar==2 & Qsortvar==2	
	qui replace FFport = 6 if mQsortvar==2 & Qsortvar==3	
	
	qui drop if missing(FFport)
	
	keep yyyymm FFport Nfirms vwret ???beta
	
	qui reshape wide vwret Nfirms ???beta, i(yyyymm) j(FFport)	
	
	qui gen double model`model'umd = (1/2) * (vwret3 + vwret6) - (1/2) * (vwret1 + vwret4)		
	rename *beta? model`model'_*beta?
	
	keep yyyymm model*
		
	* roll yyyymm forward by a month (so that yyyymm is from the viewpoint of returns, retnm)
	qui gen     yyyymm_nm = yyyymm + 1                        if mod(yyyymm,100)<12
	qui replace yyyymm_nm = 100 * (floor(yyyymm/100) + 1) + 1 if mod(yyyymm,100)==12
	qui replace yyyymm    = yyyymm_nm
	drop    yyyymm_nm				
				
	if `model'>0 {
		qui merge 1:1 yyyymm using "resmom_factors2.dta", nogenerate
		}
	save "resmom_factors2.dta", replace
	
	}
	
	

	
* Illustrate the BAB effect 
use "resmom_factors2.dta", clear

qui merge m:1 yyyymm using "fffactors.dta", keep(1 3) keepusing(mktrf smb hml rmw cma umd) nogenerate 

rename model4umd bab_capm_mktbeta
rename model8umd bab_ff5_mktbeta
replace bab_capm_mktbeta = -bab_capm_mktbeta
replace bab_ff5_mktbeta = -bab_ff5_mktbeta

regress bab_capm_mktbeta
regress bab_capm_mktbeta mktrf
regress bab_ff5_mktbeta
regress bab_ff5_mktbeta mktrf

	
	
* Compute betas for portfolios when sorted by FF5 residuals 
use "resmom_factors2.dta", clear
keep yyyymm model3*

keep if yyyymm>=197307

foreach beta in "mkt" "smb" "hml" "rmw" "cma" {
	qui gen double long_`beta'beta = (1/2) * (model3_`beta'beta3 + model3_`beta'beta6) 
	qui gen double short_`beta'beta =(1/2) * (model3_`beta'beta1 + model3_`beta'beta4)
	qui gen double ls_`beta'beta = (1/2) * (model3_`beta'beta3 + model3_`beta'beta6) - (1/2) * (model3_`beta'beta1 + model3_`beta'beta4)
	}

keep yyyymm model3umd long* short* ls*
	
foreach var of varlist ls_* {
	qui gen double se_`var' = `var'
	}
	
collapse (mean) long* short* ls* (semean) se_*

	
* Quantify the amount of residual momentum 

use "oos_tsmom_scs.dta", clear
rename tsmom? pctsmom?

merge 1:1 yyyymm using "tsmom_basic.dta", nogenerate
replace TSMom = TSMom / 100

* Merge in Fama-French factors
qui merge m:1 yyyymm using "fffactors.dta", keep(1 3) keepusing(mktrf smb hml rmw cma umd) nogenerate 

merge 1:1 yyyymm using "resmom_factors2.dta", nogenerate

rename model4umd bab_capm_mktbeta

rename model5umd bab_ff3_mktbeta
rename model6umd bab_ff3_smbbeta
rename model7umd bab_ff3_hmlbeta

rename model8umd bab_ff5_mktbeta
rename model9umd bab_ff5_smbbeta
rename model10umd bab_ff5_hmlbeta
rename model11umd bab_ff5_rmwbeta
rename model12umd bab_ff5_cmabeta

rename model13umd bab_ff5_fit

* use consistent umd (same universe)
replace umd = model0umd

drop if missing(pctsmom0)

estimates clear

forvalues model=0/3 {

	qui regress model`model'umd mktrf smb hml rmw cma
	estimates store regr`model'a
	if `model'>0 {
		qui regress model`model'umd mktrf smb hml rmw cma umd
		estimates store regr`model'b
		}
	qui regress model`model'umd mktrf smb hml rmw cma TSMom
	estimates store regr`model'c
	qui regress model`model'umd mktrf smb hml rmw cma TSMom pctsmom0
	estimates store regr`model'd		
	}

esttab regr0*, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 

esttab regr1* regr2* regr3*, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 


