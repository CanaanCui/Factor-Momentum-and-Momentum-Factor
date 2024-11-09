* Construct factors and momentum-neutral factors using the KNS signals
*
* 0. For computing daily factor returns, take just returns from the daily
*    CRSP file and add yyyymm for the *previous* month for merging
*
* 1. Open the CSV file with the characteristics, as provided by Kozak et al. 
*    (2019)
*
* 2. Shift date. Characteristics are reported as of the beginning of the month.
*    If date = 07/1963, the characteristics/weights are the 06/1963 weights.
*
* 3. Merge in (a) next month's return, (b) prior one-year return,
*    r12_2, and (c) market value of equity 
*
* 4. Impose KNS sample restrictions: remove microcaps; SCS defines microcaps as
*    follows: "To ensure that the results are not driven by very small illiquid 
*    stocks, we exclude small- cap stocks with market caps below 0.01% of 
*    aggregate stock market capitalization at each point in time."
*
* 5. Rescale weights using the Frazzini-Pedersen approach
*
* 6. Create an alternative set of weights that are orthogonal w.r.t. past return
*
* 7. Multiply weights by returns and take sums to construct the original and
*    momentum-neutral factors.
*


* REQUIRED INPUT FILES
*
* CRSPmonthly2019dec_edited.dta			Monthly CRSP file from WRDS, December 2019 vintage
* CRSPdaily2019dec_edited.dta			Daily CRSP file from WRDS, December 2019 vintage
* characteristics_anom.csv				KNS characteristics from Serhiy Kozak's website
* merged_crsp_compustat.dta				Monthly CRSP file, me = market value of equity, r12_2 = prior one year returns skipping a month, retnm = return inclusive of dividends *next month*
 
*===============================================================================
* Daily returns
*===============================================================================

* Step a. Create a file for removing stocks that are *never* eligible
use "CRSPmonthly2019dec_edited.dta", clear			

keep permno shrcd exchcd 

gen ok = 1 if shrcd>=1 & exchcd<=3 & (shrcd==10 | shrcd==11)

keep if ok==1
keep permno

duplicates drop

save "tmp_info.dta", replace

* Step b. Process daily returns

use "CRSPdaily2019dec_edited.dta", clear			

keep date permno ret

* We don't need pre-1963 returns
drop if year(date)<1963

drop if missing(ret)

* Restrict the sample to stocks that are eligible
merge m:1 permno using "tmp_info.dta", nogenerate keep(3)

* Create *lagged* yyyymm for merging
gen long yyyymm = 100 * year(date) + month(date) - 1 if month(date) > 1
replace  yyyymm = 100 * (year(date) - 1) + 12 if month(date) == 1

compress

sort permno yyyymm

save "dailyreturns.dta", replace



*===============================================================================
* 1. Open the CSV file with the KNS characteristics
*===============================================================================
import delimited using "characteristics_anom.csv", delimiter(",") clear

* Characteristics are reported as of the end of the month.
* Change the convention to end-of-month by shifting the date variable.

gen long yyyymm_nm = real(substr(date,4,4)) * 100 + real(substr(date,1,2))

gen long yyyymm = yyyymm_nm - 1 if mod(yyyymm_nm,100)>1
replace  yyyymm = (trunc(yyyymm_nm/100) - 1) * 100 + 12 if mod(yyyymm_nm,100)==1

summ yyyymm yyyymm_nm

drop date yyyymm_nm

order permno yyyymm

sort permno yyyymm


*===============================================================================
* 2. Merge in additional information 
*===============================================================================

* We want to keep the full universe of stocks to implement the market-cap filter
* for SCS replication

merge 1:1 permno yyyymm using "merged_crsp_compustat.dta", keepusing(me retnm r12_2)

keep if yyyymm>=196306

* Verify that the merge looks good

tab _merge
keep if _merge==2 | _merge==3
drop _merge


*===============================================================================
* 3. Save the file
*===============================================================================

save "KNSsignals.dta", replace


*===============================================================================
* 4. Create three sets of factors
*
*    value				value factor based on the full universe of stocks
*    value_scs 			value factor after excluding stocks with 
*                       ME<0.01% * total ME
*    value_mom			value factor with the exclusion *and* weights orthogonal
*                       to "prior one-year returns skipping a month," r12_2
*===============================================================================
local bNoDaily = 1

use "KNSsignals.dta", clear

local varctr = 0

foreach var in size value prof dur valprof fscore debtiss repurch nissa accruals growth aturnover gmargins divp ep cfp noa inv invcap igrowth sgrowth lev roaa roea sp gltnoa divg invaci mom indmom valmom valmomprof shortint mom12 momrev lrrev valuem nissm sue roe rome roa strev ivol betaarb season indrrev indrrevlv indmomrev ciss price age shvol ipo {
	
	local varctr = `varctr' + 1

	disp("Variable: `var'")
	
	use "KNSsignals.dta", clear
	
	keep permno yyyymm me `var' retnm r12_2
	
	rename `var' z
		
	*=========== Big-cap factor: drop small stocks from the definition of z2
	qui bysort yyyymm: egen totalme = sum(me)
	
	qui gen z2 = z if me >= 0.0001 * totalme & ~missing(totalme) & ~missing(me)
	qui gen z3 = z if me <  0.0001 * totalme & ~missing(totalme) & ~missing(me) 
	
	drop me totalme

	* remove stocks without characteristics
	qui drop if missing(z)
	
	* rescale z2 using the KNS/Frazzini-Pedersen formula
	qui bysort yyyymm: egen r = rank(z2)
	qui bysort yyyymm: egen n = max(r)
	qui gen rc = r / (n + 1)
	qui bysort yyyymm: egen mrc = mean(rc)
	qui gen absdev = abs(rc - mrc)
	qui bysort yyyymm: egen sumabsdev = sum(absdev)
	qui replace z2 = (rc - mrc) / sumabsdev
	drop r n rc mrc absdev sumabsdev
		
	* Compute past return-neutral weights
	* Regress raw weights on past returns; new weights are the residuals
	qui statsby _b _aR2=e(r2_a), by(yyyymm) nodots saving("regr.dta", replace): regress z2 r12_2 
	qui merge m:1 yyyymm using "regr.dta", nogenerate 
	qui gen z2m = z2 - _b_cons - _b_r12_2 * r12_2
	drop _b_cons _b_r12_2
	rename _eq2_aR2 bigR2
	
	*=========== Small-cap factor: drop big stocks from the definition of z2

	* rescale z3 using the KNS/Frazzini-Pedersen formula
	qui bysort yyyymm: egen r = rank(z3)
	qui bysort yyyymm: egen n = max(r)
	qui gen rc = r / (n + 1)
	qui bysort yyyymm: egen mrc = mean(rc)
	qui gen absdev = abs(rc - mrc)
	qui bysort yyyymm: egen sumabsdev = sum(absdev)
	qui replace z3 = (rc - mrc) / sumabsdev
	drop r n rc mrc absdev sumabsdev
		
	* Compute past return-neutral weights
	* Regress raw weights on past returns; new weights are the residuals
	qui statsby _b _aR2=e(r2_a), by(yyyymm) nodots saving("regr.dta", replace): regress z3 r12_2 
	qui merge m:1 yyyymm using "regr.dta", nogenerate 
	qui gen z3m = z3 - _b_cons - _b_r12_2 * r12_2
	drop _b_cons _b_r12_2
	rename _eq2_aR2 smallR2


	* <<MONTHLY>> Compute monthly factor returns
	preserve
	
	qui gen double z_retnm = z * retnm if ~missing(z,retnm)
	qui gen double z2_retnm = z2 * retnm if ~missing(z2,retnm)
	qui gen double z2m_retnm = z2m * retnm if ~missing(z2m,retnm)
	qui gen double z3_retnm = z3 * retnm if ~missing(z3,retnm)
	qui gen double z3m_retnm = z3m * retnm if ~missing(z3m,retnm)

	* Compute past returns (r12_2) for the stocks in the factor today
	qui gen double z2_pastret = z2 * r12_2 if ~missing(z2,z2m,r12_2)
	qui gen double z2m_pastret = z2m * r12_2 if ~missing(z2,z2m,r12_2)
	
	collapse (sum) `var'_full = z_retnm `var'_scs = z2_retnm `var'_mom = z2m_retnm `var'_scss = z3_retnm `var'_moms = z3m_retnm `var'_scspast = z2_pastret `var'_mompast = z2m_pastret (mean) `var'_bigR2=bigR2 `var'_smallR2=smallR2 (count) n1 = z_retnm n2 = z2_retnm n2m = z2m_retnm n3 = z3_retnm n3m = z3m_retnm, by(yyyymm)
	
	* Need to set return to zero if no stocks
	qui replace `var'_full = . if n1==0
	qui replace `var'_scs  = . if n2==0
	qui replace `var'_mom  = . if n2m==0
	qui replace `var'_scss = . if n3==0
	qui replace `var'_moms = . if n3m==0
	
	qui replace `var'_scspast = . if missing(`var'_mom)
	qui replace `var'_mompast = . if missing(`var'_mom)
	
	qui drop if n1==0 & n2==0 & n2m==0 & n3==0 & n3m==0
	
	rename n1 `var'_N
	rename n2 `var'_Nscs
	rename n3 `var'_Nscss
	
	drop n2m n3m
	
	* roll yyyymm forward	
	qui gen long yyyymm_nm = yyyymm + 1 if mod(yyyymm,100)<12
	qui replace  yyyymm_nm = (trunc(yyyymm/100) + 1) * 100 + 1 if mod(yyyymm,100)==12
	
	drop yyyymm
	rename yyyymm_nm yyyymm
	
	if `varctr'>1 {
		qui merge 1:1 yyyymm using "kns_factors.dta", nogenerate 
		}
	else {
		disp("new file")
		}
		
	order `var'_full `var'_scs `var'_mom `var'_scss `var'_moms `var'_scspast `var'_mompast `var'_N `var'_Nscs, last
	
	qui save "kns_factors.dta", replace
	
	restore
	
	
	* <<DAILY>> Compute daily factor returns
	if `bNoDaily'==0 {
		
		drop retnm
		
		sort permno yyyymm
		qui merge 1:m permno yyyymm using "dailyreturns.dta", nogenerate keep(1 3)
		drop yyyymm
		
		qui gen double z_ret = z * ret if ~missing(z,ret)
		qui gen double z2_ret = z2 * ret if ~missing(z2,ret)
		qui gen double z2m_ret = z2m * ret if ~missing(z2m,ret)
		qui gen double z3_ret = z3 * ret if ~missing(z3,ret)
		qui gen double z3m_ret = z3m * ret if ~missing(z3m,ret)

		* Compute past returns (r12_2) for the stocks in the factor today
		qui gen double z2_pastret = z2 * r12_2 if ~missing(z2,z2m,r12_2)
		qui gen double z2m_pastret = z2m * r12_2 if ~missing(z2,z2m,r12_2)
		
		collapse (sum) `var'_full = z_ret `var'_scs = z2_ret `var'_mom = z2m_ret `var'_scss = z3_ret `var'_moms = z3m_ret `var'_scspast = z2_pastret `var'_mompast = z2m_pastret (count) n1 = z_ret n2 = z2_ret n2m = z2m_ret n3 = z3_ret n3m = z3m_ret, by(date)
		
		* Need to set return to zero if no stocks
		qui replace `var'_full = . if n1==0
		qui replace `var'_scs  = . if n2==0
		qui replace `var'_mom  = . if n2m==0
		qui replace `var'_scss = . if n3==0
		qui replace `var'_moms = . if n3m==0
		
		qui replace `var'_scspast = . if missing(`var'_mom)
		qui replace `var'_mompast = . if missing(`var'_mom)
		
		qui drop if n1==0 & n2==0 & n2m==0 & n3==0 & n3m==0
		
		rename n1 `var'_N
		rename n2 `var'_Nscs
		rename n3 `var'_Nscss
		
		drop n2m n3m
			
		if `varctr'>1 {
			qui merge 1:1 date using "kns_factors_daily.dta", nogenerate 
		}
		else {
			disp("new file")
		}
			
		order `var'_full `var'_scs `var'_mom `var'_scss `var'_moms `var'_scspast `var'_mompast `var'_N `var'_Nscs, last
		
		qui save "kns_factors_daily.dta", replace
		
	}
	
}
	
	
*===============================================================================
* Verify that the replication of the SCS factors is good
*
* Managed_portfolios_anom_50.csv is the file with factor returns provided by 
* Kozak et al. on Serhiy Kozak's website
*
* The correlation between r_value and value_scs should be high--value_scs uses 
* the same universe
*===============================================================================

import delimited using "managed_portfolios_anom_50.csv", delimiter(",") clear

gen long yyyymm = real(substr(date,4,4)) * 100 + real(substr(date,1,2))

order yyyymm

merge 1:1 yyyymm using "kns_factors.dta", nogenerate

order r_size size_scs size_mom r_value value_scs value_mom r_prof prof_scs prof_mom 

corr r_size size_scs
corr r_size size_mom
corr r_value value_scs
corr r_value value_mom
corr r_prof prof_scs
corr r_prof prof_mom

*===============================================================================
* 1. Count the predictors
* 2. Illustrate that the factors are past-return neutral
*===============================================================================

use "kns_factors.dta", clear

local N=0
foreach var of varlist *_scs {
	local N = `N'+1
	disp("`var'")
}
disp("Final count: `N'")

keep yyyymm value_scs value_mom value_scspast value_mompast betaarb_scs betaarb_mom betaarb_scspast betaarb_mompast

gen avalue_scspast = abs(value_scspast)
gen avalue_mompast = abs(value_mompast)
gen abetaarb_scspast = abs(betaarb_scspast)
gen abetaarb_mompast = abs(betaarb_mompast)

summ avalue* abetaarb* if yyyymm<=201901


*===============================================================================
* How much of the variation is explained by past returns?
*===============================================================================

use "kns_factors.dta", clear

keep yyyymm *_bigR2
local factorid = 0

foreach var in size value prof dur valprof fscore debtiss repurch nissa accruals growth aturnover gmargins divp ep cfp noa inv invcap igrowth sgrowth lev roaa roea sp gltnoa divg invaci shortint lrrev valuem nissm sue roe rome roa strev ivol betaarb season indrrev indrrevlv ciss price age shvol ipo {
	
	local factorid = `factorid' + 1
	rename `var'_bigR2 r2factor`factorid'
	gen factorname`factorid' = "`var'"
	
	}

keep yyyymm r2factor* factorname* 

reshape long r2factor factorname, i(yyyymm) j(factorid)

collapse (mean) r2factor, by(factorid factorname)


*===============================================================================
* Extract principal components and measure the amount of time-series factor
* momentum
*
* Compute time-series factor momentum returns out-of-sample
*
* Initial step: Create simple monthly and daily files with just the factors
*
* Step 1. Estimate PCs from daily data up to month t
* Step 2. Compute scores using monthly data up to month t+1
*         Scale PC variances (up to month t) to match the average variance
*         of raw factors (up to month t)
* Step 3. Compute tsmom return for month t+1
*===============================================================================

* CHOOSE FACTOR STYLE:
* 
* 1) scs = standard factors
* 2) mom = momentum-neutral factors

local factor_style = "mom"


* Prepare factors: Monthly
local factorctr = 0

use "kns_factors.dta", clear

foreach var in size value prof dur valprof fscore debtiss repurch nissa accruals growth aturnover gmargins divp ep cfp noa inv invcap igrowth sgrowth lev roaa roea sp gltnoa divg invaci shortint lrrev valuem nissm sue roe rome roa strev ivol betaarb season indrrev indrrevlv ciss price age shvol ipo {
		
	local factorctr = `factorctr'+1
		
	rename `var'_`factor_style' factor`factorctr'
		
}

keep yyyymm factor* 

* replace missing factor returns with zeros - this seems to be implied by footnote 16 in the KNS2020

qui mvencode _all, mv(0) override

save "kns_factors_clean.dta", replace


* Prepare factors: Daily
local factorctr = 0

use "kns_factors_daily.dta", clear

foreach var in size value prof dur valprof fscore debtiss repurch nissa accruals growth aturnover gmargins divp ep cfp noa inv invcap igrowth sgrowth lev roaa roea sp gltnoa divg invaci shortint lrrev valuem nissm sue roe rome roa strev ivol betaarb season indrrev indrrevlv ciss price age shvol ipo {
		
	local factorctr = `factorctr'+1
		
	rename `var'_`factor_style' factor`factorctr'
		
}

keep date factor* 

qui mvencode _all, mv(0) override

save "kns_factors_daily_clean.dta", replace




local loop_ctr = 0

forvalues yyyy=1973/2020 {
	forvalues mo=1/12 {
		
		local loop_ctr = `loop_ctr' + 1
		
		* Step 1: Compute PCs from daily data
		* ======
		
		use "kns_factors_daily_clean.dta", clear
		
		* Last day of the month?
		if `mo'<12 {
			local lastday = mdy(`mo'+1,1,`yyyy') - 1
		}
		else {
			local lastday = mdy(1,1,`yyyy'+1) - 1
			}
		
		qui keep if date<=`lastday'

		* VERIFY
		*sort date
		*list date in -1 
		
		qui pca factor*
		
		* Step 2: Compute scores using monthly data
		* ======
				
		use "kns_factors_clean.dta", clear
		
		if `mo'<12 {
			local oos_month = 100 * `yyyy' + `mo' + 1
		}
		else {
			local oos_month = 100 * (`yyyy' + 1) + 1
		}

		qui keep if yyyymm<=`oos_month'
		
		disp("`oos_month'")
		
		* compute average variance of factors
		preserve
		
		* remove oos month
		qui drop if yyyymm==`oos_month'
		
		collapse (sd) factor*
		foreach var of varlist factor* {
			qui replace `var' = `var'^2
			}
		egen avg = rowmean(factor*)
		qui summ avg
		local avgvar = `r(mean)'

		restore

		* compute scores

		qui predict pc1-pc47

		keep yyyymm pc*

		* rescale to match average variance of the actual factors
		foreach var of varlist pc* {
			qui summ `var'	
			qui replace `var' = `var' * sqrt(`avgvar') / `r(sd)'
		}
		
		
		* Step 3: Construct factor momentum strategies
		* ======
		
		qui reshape long pc, i(yyyymm) j(factorid)

		gen long monthindex = (trunc(yyyymm/100)-1900)*12 + mod(yyyymm,100)

		qui xtset factorid monthindex

		* compute returns over the prior year
		qui gen ret12 = L1.pc + L2.pc + L3.pc + L4.pc + L5.pc + L6.pc + L7.pc + L8.pc + L9.pc + L10.pc + L11.pc + L12.pc

		* we only need the one OOS month
		qui keep if yyyymm==`oos_month'
		
		qui keep if ~missing(ret12)

		gen tsmom = sign(ret12) * pc	
		
		gen byte group = trunc((factorid - 1) / 10) + 1
		
		* compute group-specific strategies
		preserve

		collapse (mean) tsmom, by(group yyyymm)
				
		qui reshape wide tsmom, i(yyyymm) j(group)
		
		qui save "oos_tsmom_tmp.dta", replace
		restore
		
		* all
		collapse (mean) tsmom0 = tsmom, by(yyyymm)
		
		qui merge 1:1 yyyymm using "oos_tsmom_tmp.dta", nogenerate
		
		if `loop_ctr'>1 {
			append using "oos_tsmom_`factor_style'.dta"
		}
		qui save "oos_tsmom_`factor_style'.dta", replace
			
	}
}


	
	
	
	
*===============================================================================
* Compute time-series factor momentum for KNS factors and momentum-neutral KNS 
* factors
*
* Do not include factors related to momentum
*===============================================================================

use "kns_factors.dta", clear

local factorctr = 0

foreach var in size value prof dur valprof fscore debtiss repurch nissa accruals growth aturnover gmargins divp ep cfp noa inv invcap igrowth sgrowth lev roaa roea sp gltnoa divg invaci shortint lrrev valuem nissm sue roe rome roa strev ivol betaarb season indrrev indrrevlv ciss price age shvol ipo {
	
	local factorctr = `factorctr' + 1

	drop `var'_full 
	
	rename `var'_scs factor`factorctr'
	rename `var'_mom mom_factor`factorctr'
	
}
	
keep yyyymm factor* mom_factor*
	
save "temp_factors.dta", replace


forvalues yyyy=1973/2020 {
	forvalues mo=1/12 {
		
		local loop_ctr = `loop_ctr' + 1
		
		if `mo'<12 {
			local oos_month = 100 * `yyyy' + `mo' + 1
		}
		else {
			local oos_month = 100 * (`yyyy' + 1) + 1
		}
		
		use "temp_factors.dta", clear
		
		* compute average variances of factors
		qui drop if yyyymm>`oos_month'

		disp("`oos_month'")
		
		preserve
		
		qui drop if yyyymm==`oos_month'
		
		collapse (sd) factor* mom_factor*
		foreach var of varlist factor* mom_factor* {
			qui replace `var' = `var'^2
		}
		egen avg = rowmean(factor*)
		egen mom_avg = rowmean(mom_factor*)
		qui summ avg
		local avgvar = `r(mean)'
		qui summ mom_avg
		local mom_avgvar = `r(mean)'
		
		restore
		
		* rescale to match average variance of the actual factors
		foreach var of varlist factor* {
			qui summ `var'	
			if `r(N)'>12 {
				qui replace `var' = `var' * sqrt(`avgvar') / `r(sd)'
				}
			}

		foreach var of varlist mom_factor* {
			qui summ `var'	
			if `r(N)'>12 {
				qui replace `var' = `var' * sqrt(`mom_avgvar') / `r(sd)'
			}
		}


		qui reshape long factor mom_factor, i(yyyymm) j(factorid)
			
		gen long dateindex = (trunc(yyyymm/100) - 1900) * 12 + mod(yyyymm,100)
		qui xtset factorid dateindex

		* Compute prior one-year returns for the factors
		foreach var of varlist factor mom_factor {
			qui gen L`var' = L1.`var' + L2.`var' + L3.`var' + L4.`var' + L5.`var' + L6.`var' + L7.`var' + L8.`var' + L9.`var' + L10.`var' + L11.`var' + L12.`var'
		}

		* we only need the one OOS month
		qui keep if yyyymm==`oos_month'
			
		gen     tsmom = sign(Lfactor) * factor 
		gen mom_tsmom = sign(Lmom_factor) * mom_factor
		
		collapse (mean) tsmom mom_tsmom, by(yyyymm)	

		if `loop_ctr'>1 {
			append using "oos_tsmom_raw.dta"
		}
		qui save "oos_tsmom_raw.dta", replace		
			
	}
}		

			
			


*===============================================================================
*
* Table 7: Factor momentum in high- and low-eigenvalue factors
* ========
*
* 1. Amount of momentum in raw factors
* 2. Amount of momentum in pc-based factors
* 3. Do pc-based factors span raw factors, and vice versa?
*
*===============================================================================

use "oos_tsmom_raw.dta", clear
rename tsmom tsmom_raw 
rename mom_tsmom mom_tsmom_raw

merge 1:1 yyyymm using "oos_tsmom_mom.dta", nogenerate
rename tsmom? pcmom_tsmom?

merge 1:1 yyyymm using "oos_tsmom_scs.dta", nogenerate
rename tsmom? pctsmom?

* Merge in Fama-French factors
qui merge m:1 yyyymm using "fffactors.dta", keep(1 3) keepusing(mktrf smb hml rmw cma umd) nogenerate 

foreach var of varlist tsmom_raw-umd {
	qui replace `var' = 100 * `var'
}
	
keep if yyyymm>=197307

estimates clear

qui regress tsmom_raw mktrf smb hml rmw cma
estimates store regr1
qui regress tsmom_raw mktrf pctsmom0 smb hml rmw cma
estimates store regr2
qui regress tsmom_raw  umd mktrf  smb hml rmw cma
estimates store regr3

qui regress pctsmom0 mktrf smb hml rmw cma
estimates store regr4
qui regress pctsmom0 mktrf tsmom_raw smb hml rmw cma
estimates store regr5
qui regress pctsmom0 umd mktrf  smb hml rmw cma
estimates store regr6

forvalues group=1/5 {
	
	qui regress pctsmom`group' mktrf smb hml rmw cma
	estimates store regra`group'

}

esttab regr? regra?, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 



	
*===============================================================================
*
* Table 8: Factor momentum in momentum-neutral factors
*
*===============================================================================

estimates clear

qui regress tsmom_raw mktrf smb hml rmw cma
estimates store regr1
qui regress tsmom_raw mktrf smb hml rmw cma mom_tsmom_raw
estimates store regr2
qui regress mom_tsmom_raw mktrf smb hml rmw cma
estimates store regr3
qui regress mom_tsmom_raw mktrf smb hml rmw cma tsmom_raw
estimates store regr4
qui regress pctsmom0 mktrf smb hml rmw cma
estimates store regr5
qui regress pctsmom0 mktrf smb hml rmw cma pcmom_tsmom0 
estimates store regr6
qui regress pcmom_tsmom0 mktrf smb hml rmw cma
estimates store regr7
qui regress pcmom_tsmom0 mktrf smb hml rmw cma pctsmom0 
estimates store regr8

esttab regr?, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 


estimates clear

forvalues group=0/5 {
	
	qui regress pcmom_tsmom`group' mktrf
	estimates store regra`group'

	qui regress pcmom_tsmom`group' mktrf umd
	estimates store regrb`group'

	qui regress umd pcmom_tsmom`group' mktrf 
	estimates store regrc`group'
	
}

esttab regra? regrb?, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 

esttab regrc?, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 


	
	
*===============================================================================
*
* FMOM2's t-values for Figure 1
*
* FMOM2 in FF5
* FMOM2 in FF5 + 5 UMDs
*
*===============================================================================

use "otherumds.dta", clear
gen long yyyymm = 100 * year + month

merge 1:1 yyyymm using "oos_tsmom_scs.dta", nogenerate keepusing(tsmom0)
rename tsmom0 FMOM2

* Merge in Fama-French factors
qui merge m:1 yyyymm using "fffactors.dta", keep(1 3) keepusing(mktrf smb hml rmw cma umd) nogenerate 

keep if yyyymm>=196307 & yyyymm<=201912

summ yyyymm if ~missing(FMOM2)

estimates clear

qui regress FMOM2 mktrf smb hml rmw cma
estimates store regr1

qui regress FMOM2 mktrf smb hml rmw cma umd UMD_Industry UMD_Sharpe UMD_Novy UMD_IndAdj
estimates store regr2

esttab regr?, b(5) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar 


	
	
