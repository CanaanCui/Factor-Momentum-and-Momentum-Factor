*===============================================================================
* Calculation of other forms of momentum
*===============================================================================

/**************************Industry neutral Momentum**************************/
use "CRSPMonthly1926-2019.dta", clear

gen me=abs(prc* shrout)
gen time=ym(year,month)
xtset permno time

* momentum is the average of t-12 to t-2 returns
tssmooth ma MA=ret ,window(11)
gen MOM=l.MA

* Fama-French 49 industry-neutral
ffind siccd, newvar(ff49code) type(49)
keep if !missing(ff49code)
gen sortvar=MOM

* Sort
qui bysort time ff49code: egen m  = mean(sortvar) if ~missing(ff49code)
qui gen double sortvar_ff = sortvar - m

xtset
gen lag_me=l.me
drop if sortvar_ff==. | lag_me==. | ret==.

qui gen double nyse_sortvar_ff = sortvar_ff if exchcd==1 
qui bysort time: egen sv30 = pctile(nyse_sortvar_ff), p(30)
qui bysort time: egen sv70 = pctile(nyse_sortvar_ff), p(70)

qui gen byte Qsortvar_ff = 1 if                     sortvar_ff <= sv30 & ~missing(sortvar_ff)
qui replace  Qsortvar_ff = 2 if sortvar_ff > sv30 & sortvar_ff <= sv70 & ~missing(sortvar_ff)
qui replace  Qsortvar_ff = 3 if sortvar_ff > sv70                      & ~missing(sortvar_ff)

qui gen double nyse_size = me if exchcd==1 
qui bysort time: egen size50 = pctile(nyse_size), p(50)

qui gen byte Qsortsize = 1 if me <= size50 & ~missing(me)
qui replace  Qsortsize = 2 if me > size50  & ~missing(me)
		
* 6 value weighted portfolios at the intersection of 3 momentum and 2 size
egen IndAdjMom= wtmean(ret), by(Qsortvar_ff Qsortsize time) weight(lag_me)

collapse (mean) IndAdjMom ,by(Qsortvar_ff Qsortsize year month)

save "Temp/Mom_IndustryNeutral.dta", replace

/************************** NovyMarx -12 to -7 Momentum**************************/
use "CRSPMonthly1926-2019.dta", clear

gen me=abs(prc* shrout)
gen time=ym(year,month)
xtset permno time

* momentum is the average of t-12 to t-7 returns
tssmooth ma MA=ret ,window(6)
gen MOM=l6.MA	

gen sortvar=MOM

xtset
gen lag_me=l.me
drop if sortvar==. | lag_me==. | ret==.

qui gen double nyse_sortvar= sortvar if exchcd==1 
qui bysort time: egen sv30 = pctile(nyse_sortvar), p(30)
qui bysort time: egen sv70 = pctile(nyse_sortvar), p(70)

qui gen byte Qsortvar= 1 if                     sortvar <= sv30 & ~missing(sortvar)
qui replace  Qsortvar= 2 if sortvar> sv30 & sortvar  <= sv70 & ~missing(sortvar)
qui replace  Qsortvar= 3 if sortvar> sv70                      & ~missing(sortvar)

qui gen double nyse_size = me if exchcd==1 
qui bysort time: egen size50 = pctile(nyse_size), p(50)

qui gen byte Qsortsize = 1 if me <= size50 & ~missing(me)
qui replace  Qsortsize = 2 if me > size50  & ~missing(me)
		
* 6 value weighted portfolios at the intersection of 3 momentum and 2 size
egen NovyMom= wtmean(ret), by(Qsortvar Qsortsize time) weight(lag_me)

collapse (mean) NovyMom ,by(Qsortvar Qsortsize year month)
save "Temp/Mom_Novy.dta", replace

/************************** Sharpe ratio Momentum **************************/	
use "CRSPMonthly1926-2019.dta", clear

gen me=abs(prc* shrout)
gen time=ym(year,month)
xtset permno time

tssmooth ma MA=ret ,window(11)
gen MOM=l.MA	

/* compute t-12 to t-2 standard deviation of residual returns*/
tsfill
mvsumm ret, stat(sd) win(11) gen(SDt11_t1) force end
gen Sharpe_t12t2=MOM/l.SDt11_t1
drop if Sharpe_t12t2==.

gen sortvar=Sharpe_t12t2

xtset
gen lag_me=l.me
drop if sortvar==. | lag_me==. | ret==.

qui gen double nyse_sortvar= sortvar if exchcd==1 
qui bysort time: egen sv30 = pctile(nyse_sortvar), p(30)
qui bysort time: egen sv70 = pctile(nyse_sortvar), p(70)

qui gen byte Qsortvar= 1 if                     sortvar <= sv30 & ~missing(sortvar)
qui replace  Qsortvar= 2 if sortvar> sv30 & sortvar  <= sv70 & ~missing(sortvar)
qui replace  Qsortvar= 3 if sortvar> sv70                      & ~missing(sortvar)

qui gen double nyse_size = me if exchcd==1 
qui bysort time: egen size50 = pctile(nyse_size), p(50)

qui gen byte Qsortsize = 1 if me <= size50 & ~missing(me)
qui replace  Qsortsize = 2 if me > size50  & ~missing(me)

* 6 value weighted portfolios at the intersection of 3 momentum and 2 size
egen Sharpe= wtmean(ret), by(Qsortvar Qsortsize time) weight(lag_me)

collapse (mean) Sharpe ,by(Qsortvar Qsortsize year month)
save "Temp/Mom_Sharpe.dta", replace

/**************************Grinblatt & Moskowitz 1999 Industry Momentum**************************/

use "CRSPMonthly1926-2019.dta", clear
gen me=abs(prc* shrout)
gen time=ym(year,month)
xtset permno time

* Fama-French 17 industry-neutral
ffind siccd, newvar(ff17code) type(17)
xtset
gen lag_me=l.me
drop if ff17code==. | lag_me==. | ret==.
egen IndRet= wtmean(ret), by(ff17code time) weight(lag_me)
collapse (mean) IndRet, by(ff17code year month time)

xtset ff17code time
tssmooth ma MOM=IndRet,window(12)

* buy top 3 and short top 3 
sort time MOM
qui bysort time: gen Rank= _n

qui gen Qsortvar = 1 if (Rank==1 | Rank==2 | Rank==3)
qui replace Qsortvar = 2 if (Rank==15 | Rank==16 | Rank==17)

collapse (mean) IndRet if !missing(Qsortvar),by(Qsortvar year month)
save "Temp/Mom_IndustryRet.dta", replace

*===============================================================================
* merge momentums in one file	
	
* Industry-adjusted momentum
use "Temp/Mom_IndustryNeutral.dta", clear
collapse (mean) IndAdjMom ,by(Qsortvar_ff year month)
reshape wide IndAdjMom, i(year month) j(Qsortvar_ff)
gen UMD_IndAdj=100*(IndAdjMom3-IndAdjMom1)
keep year month UMD_IndAdj
save "otherumds.dta", replace

* Intermediate momentum of Novy-Marx
use "Temp/Mom_Novy.dta", clear
collapse (mean) NovyMom ,by(Qsortvar year month)
reshape wide NovyMom ,i(year month) j(Qsortvar)
gen UMD_Novy=100*(NovyMom3-NovyMom1)
keep year month UMD_Novy
merge 1:1 year month using "otherumds.dta", nogenerate
save "otherumds.dta", replace

* Sharpe ratio momentum
use "Temp/Mom_Sharpe.dta", clear
collapse (mean) Sharpe ,by(Qsortvar year month)
reshape wide Sharpe ,i(year month) j(Qsortvar)
gen UMD_Sharpe=100*(Sharpe3- Sharpe1)
keep year month UMD_Sharpe
merge 1:1 year month using "otherumds.dta", nogenerate
save "otherumds.dta", replace

* Grinblatt Moskiwitz Indsutry momentum
use "Temp/Mom_IndustryRet.dta", clear
reshape wide IndRet ,i(year month) j(Qsortvar)
gen UMD_Industry=100*(IndRet2-IndRet1)
keep year month UMD_Industry
merge 1:1 year month using Tests/NewDataSets/otherumds.dta, nogenerate
save "otherumds.dta", replace
