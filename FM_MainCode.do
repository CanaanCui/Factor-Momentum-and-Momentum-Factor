*===============================================================================
* Choose the data folder:
*===============================================================================

cd "/Users/canaandrinkwater/Desktop/Thesis 3/References/fmmf/jofi13131-sup-0002-replicationcode/Data"

*===============================================================================
* Table I: Means, t-values, and standard deviation for factors 
*===============================================================================
use "anomalies.dta", clear

gen long yyyymm = 100 * year + month

collapse (mean) ret (sd) sd=ret (semean) se=ret (count) n=ret (min) yyyymm, by(global anomaly)	
	gen tstat=ret/se
	replace ret=ret*12
	replace sd=sd*sqrt(12)
	
format ret sd tstat %12.2f 
keep anomaly yyyymm ret sd tstat 
order anomaly yyyymm ret sd tstat 

*===============================================================================
* Table II and AI: factor returns conditional of their past returns 
*===============================================================================
use "anomalies.dta", clear

egen Pan=group(global anomaly)
xtset Pan time
by Pan,sort:gen n=_n

tssmooth ma MA=ret ,window(12)
drop if n<13

gen flag=sign(MA)
replace flag=0 if flag==-1

* Individual time series regressions of factor returns on their past 12 month returns
foreach v in "a0" "a0t" "b0" "b0t" "R20" "a1" "a1t" "b1" "b1t" "R21"  {
	qui gen `v' = .
	}
	
* Duplicate the data set so that we can run pooled regressions as Pan==23
preserve
qui sum Pan
replace Pan=`r(max)'+1
replace anomaly="Everything"
save "Temp/tmp.dta", replace
restore
append using "Temp/tmp.dta"

qui sum Pan
forvalues i=1/`r(max)' {

	* Spec 1: x = average return
	if `i'<`r(max)'+1 {
		qui reg ret  MA if Pan==`i'
		}
	else {
		qui reg ret  MA if Pan==`i', cluster(time)
		}
	
	qui replace a0 = _b[_cons] if Pan==`i'
	qui replace a0t = _b[_cons] / _se[_cons] if Pan==`i'
	qui replace b0 = _b[MA] if Pan==`i'
	qui replace b0t = _b[MA] / _se[MA] if Pan==`i'

	* Spec 2: x = average return > 0
	if `i'<`r(max)'+1 {
		qui reg ret flag if Pan==`i'
		}
	else {
		qui reg ret flag if Pan==`i', cluster(time)
		}
	
	qui replace a1 = _b[_cons] if Pan==`i'
	qui replace a1t = _b[_cons] / _se[_cons] if Pan==`i'
	qui replace b1 = _b[flag] if Pan==`i'
	qui replace b1t = _b[flag] / _se[flag] if Pan==`i'

	}

collapse (mean) a0 a0t b0 b0t a1 a1t b1 b1t , by(global Pan anomaly)
  
	
*===============================================================================
* Table III. Factor Momentum in High- and Low-Eigenvalue Factors
*===============================================================================
* Panel A - total sample statistics 	
use "oos_tsmom_scs.dta",clear

foreach var of varlist tsmom* {
	replace `var'=`var'*100
	gen sd_`var'=`var'
}

collapse (mean) tsmom* (sd) sd* (count) N=tsmom1
reshape long tsmom sd_tsmom,i(N) j(subset)
gen sharpe=tsmom/sd_tsmom
gen tstat=sqrt(N)*tsmom/sd_tsmom
list subset tsmom tstat

* Panel A - statistics in the first half and second half
use "oos_tsmom_scs.dta",clear

sum yyyymm,d
gen period=1 if yyyymm<r(p50)
replace period=2 if yyyymm>=r(p50)

foreach var of varlist tsmom* {
	replace `var'=`var'*100
	gen sd_`var'=`var'
}

collapse (mean) tsmom* (sd) sd* (count) N=tsmom1,by(period)
reshape long tsmom sd_tsmom,i(period) j(subset)
gen sharpe=tsmom/sd_tsmom
gen tstat=sqrt(N)*tsmom/sd_tsmom
list subset period tsmom tstat,sep(0)
	
* Panel B and C - spanning tests
use "oos_tsmom_scs.dta",clear
merge 1:1 yyyymm using "fffactors.dta", nogenerate keep(3)
sum yyyymm,d
gen period=1 if yyyymm<r(p50)
replace period=2 if yyyymm>=r(p50)
gen x=1

estimates clear
forvalues i=2/5 {
	qui reg tsmom1 mktrf smb hml rmw cma tsmom`i' period#x,nocons
	qui estimates store regmain`i'
	qui reg tsmom`i' mktrf smb hml rmw cma tsmom1 period#x,nocons
	qui estimates store regother`i'
}

qui reg tsmom1 mktrf smb hml rmw cma tsmom2- tsmom5 period#x,nocons
qui estimates store regmain6

* Panel B: Explaining factor momentum in low-eigenvalue PC factors
esttab regother*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar
	
* Panel C: Explaining factor momentum in high-eigenvalue PC factors
esttab regmain*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar


*===============================================================================
* Create a simple TS factor based on the off-the-shelf factors
*===============================================================================

use "anomalies.dta", clear
drop if (ano=="umd" | ano=="glumd" )
gen yyyymm=year*100+month

egen Pan=group(global anomaly)
xtset Pan time
by Pan,sort:gen n=_n

set more off

tssmooth ma MA=ret ,window(12)
drop if n<13

gen sign=sign(MA)

/* TSMOM */
gen tsmomret=sign*ret
by time,sort:egen TSMom=mean(tsmomret)
		
collapse (mean) TSMom ,by(year month yyyymm)

save "Temp/TSFactor.dta", replace


*===============================================================================
* Table IV. Pricing Momentum-Sorted Portfolios with Momentum and Factor Momentum
*===============================================================================
use "P10UMD.dta", clear

merge 1:1 year month using "FactorFF5.dta",keep(3) nogenerate
merge 1:1 year month using "Temp/TSFactor.dta",keep(3) nogenerate
merge 1:1 year month using "FactorUMD.dta",keep(3) nogenerate

merge 1:1 yyyymm using "oos_tsmom_scs.dta", nogenerate
rename tsmom* pctsmom*

* excess returns 
forvalues i=1/10{
	gen ExcessP`i'=p`i'-rf
}

* portfolio 11 is high-minus-low
gen ExcessP11=ExcessP10-ExcessP1

estimates clear

* explain 10 momentum sorted portfolio with FF5
forvalues i=1/11 {
	qui reg ExcessP`i' mktrf smb hml cma rmw
	estimates store regr`i'
}
esttab regr*, b(3) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain 10 momentum sorted portfolio with FF5+UMD
estimates clear
forvalues i=1/11 {
	qui reg ExcessP`i' mktrf smb hml cma rmw umd
	estimates store regr`i'
}
esttab regr*, b(3) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain 10 momentum sorted portfolio with FF5+TSMOM
estimates clear
forvalues i=1/11 {
	qui reg ExcessP`i' mktrf smb hml cma rmw TSMom
	estimates store regr`i'
}
esttab regr*, b(3) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain 10 momentum sorted portfolio with FF5+PCMOM
replace pctsmom1=pctsmom1*100
estimates clear
forvalues i=1/11 {
	qui reg ExcessP`i' mktrf smb hml cma rmw pctsmom1
	estimates store regr`i'
}

esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* GRS tests 
drop ExcessP11

grstest2 ExcessP*,flist(mktrf  smb  hml cma rmw)
grstest2 ExcessP*,flist(mktrf  smb  hml cma rmw umd)
grstest2 ExcessP*,flist(mktrf  smb  hml cma rmw TSMom)
grstest2 ExcessP*,flist(mktrf  smb  hml cma rmw pctsmom1)

*===============================================================================
* Table V. Alternative Definitions of Momentum: Spanning Tests
*===============================================================================

use "otherumds.dta", clear
merge 1:1 year month using "FactorUMD.dta",keep(3) nogenerate
merge 1:1 year month using "FactorFF5.dta",keep(3) nogenerate
merge 1:1 year month using "Temp/TSFactor.dta",keep(3) nogenerate

merge 1:1 yyyymm using "oos_tsmom_scs.dta", nogenerate

rename tsmom* pctsmom*
replace pctsmom1=pctsmom1*100

* summary statistics - Panel A 

estimates clear
local j=0

foreach momvar of varlist umd UMD_IndAdj UMD_Industry UMD_Novy UMD_Sharpe TSMom pctsmom1 {
	
	local j=`j'+1

	qui reg `momvar' mktrf smb hml rmw cma
	estimates store regr`j'
		
}

esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain other momentums with off the shelf factor momentum
estimates clear
foreach momvar of varlist umd UMD_IndAdj UMD_Industry UMD_Novy UMD_Sharpe   {
	
	local j=`j'+1

	qui reg `momvar' mktrf smb hml rmw cma TSMom
	estimates store regr`j'
		
}

esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain other momentums with PC1-PC10

estimates clear
foreach momvar of varlist umd UMD_IndAdj UMD_Industry UMD_Novy UMD_Sharpe   {
	
	local j=`j'+1

	qui reg `momvar' mktrf smb hml rmw cma pctsmom1
	estimates store regr`j'
		
}

esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain off the shelf factor momentum using other momentums 
estimates clear
foreach momvar of varlist umd UMD_IndAdj UMD_Industry UMD_Novy UMD_Sharpe   {
	
	local j=`j'+1

	qui reg TSMom `momvar' mktrf smb hml rmw cma 
	estimates store regr`j'
		
}

qui reg TSMom umd UMD* mktrf smb hml rmw cma 
estimates store regr6
	
esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar

* explain PC1-PC10 momentum using other momentums 

estimates clear
foreach momvar of varlist umd UMD_IndAdj UMD_Industry UMD_Novy UMD_Sharpe   {
	
	local j=`j'+1

	qui reg pctsmom1 `momvar' mktrf smb hml rmw cma 
	estimates store regr`j'
		
}

qui reg pctsmom1 umd UMD* mktrf smb hml rmw cma 
estimates store regr6
	
esttab regr*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar



*===============================================================================
* Table VII. Residual Momentum versus Factor Momentum: Actual Data
*===============================================================================
use "oos_tsmom_scs.dta",clear

rename tsmom* pctsmom*

merge 1:1 yyyymm using "Temp/TSFactor.dta",keep(3) nogenerate
replace TSMom = TSMom / 100

* Merge in Fama-French factors
qui merge m:1 yyyymm using "FactorFF5.dta", keep(1 3) keepusing(mktrf smb hml rmw cma) nogenerate 

merge 1:1 yyyymm using "resmom_factors2.dta", nogenerate

* use consistent umd (same universe)
gen umd = model0umd

rename model8umd bab_ff5_mktbeta
rename model9umd bab_ff5_smbbeta
rename model10umd bab_ff5_hmlbeta
rename model11umd bab_ff5_rmwbeta
rename model12umd bab_ff5_cmabeta

drop if missing(model1umd) | missing(TSMom) 

* univariate spanning tests of residual momentum 

foreach var of varlist pctsmom* TSMom mktrf smb hml rmw cma umd bab* model* {
	replace `var'=`var'*100
}

estimates clear

forvalues model=0/3 {
	qui regress model`model'umd 
	estimates store regr`model'1

	qui regress model`model'umd TSMom
	estimates store regr`model'2

	qui regress model`model'umd pctsmom1
	estimates store regr`model'3
}

forvalues model=0/3 {
	esttab regr`model'*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar
}

* Internet appendix: Table IA.VIII. Residual Momentum Versus Factor Momentum: Alternative Factor Models
* multivariate tests of residual momentum 
order bab_ff5_mktbeta bab_ff5_smbbeta bab_ff5_hmlbeta bab_ff5_rmwbeta bab_ff5_cmabeta
keep if !missing(pctsmom1)
estimates clear

forvalues model=0/3 {
* 1: Control is FF5 
	qui regress model`model'umd 
	estimates store regr`model'0
	
	qui regress model`model'umd mktrf smb hml rmw cma
	estimates store regr`model'1

	qui regress model`model'umd mktrf smb hml rmw cma pctsmom1
	estimates store regr`model'2

	qui regress model`model'umd mktrf smb hml rmw cma bab* pctsmom1
	estimates store regr`model'4
	
}

forvalues model=0/3 {
	esttab regr`model'*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar
}


*===============================================================================
* Table VIII. Unconditional and Conditional Correlations with the Momentum Factor
*===============================================================================
use "anomalies.dta", clear

quietly drop if (ano=="umd" | ano=="glumd" )

qui egen Pan=group(global anomaly)
qui xtset Pan time
qui set more off
qui by Pan,sort:gen n=_n
qui xtset
tssmooth ma MA=ret ,window(12)
drop if n<13
gen CondRet=MA*ret 
rename ret AnomalyRet
merge m:1 year month using "FactorUMD.dta", nogenerate keep(3)

gen mypan=.
gen rawcorr=.
gen rawN=.
gen upcorr=.
gen upN=.
gen downcorr=.
gen downN=.
gen condcorr=.
gen condN=.

forvalues i=1/20{
	qui replace mypan=`i' in `i'
	qui corr AnomalyRet umd if Pan==`i'
	qui replace rawcorr=r(rho) in `i'
	qui replace rawN=r(N) in `i'
	qui corr AnomalyRet umd if Pan==`i' & MA>0
	qui replace upcorr=r(rho) in `i'
	qui replace upN=r(N) in `i'
	qui corr AnomalyRet umd if Pan==`i' & MA<0
	qui replace downcorr=r(rho) in `i'
	qui replace downN=r(N) in `i'
	qui corr CondRet umd if Pan==`i'
	qui replace condcorr=r(rho) in `i'
	qui replace condN=r(N) in `i'
}
		
/* Fisher Z-test, H0: correlations are equal */
gen Z=.
gen pvalue=.

forvalues i=1/20 {
	scalar mu_Z = atanh( condcorr[`i'] ) - atanh( rawcorr[`i'] )
	scalar sigma_Z = sqrt(1/( rawN[`i'] -3)+1/( rawN[`i'] -3))
	qui replace Z= mu_Z/sigma_Z in `i'
	qui replace pvalue= 2*normal(-abs(Z)) in `i'
}
		
	/* unbalanced Fisher Z-test for up and down, H0: correlations are equal*/
	gen Z_un=.
	gen pvalue_un=.

	forvalues i=1/20{
		scalar mu_Z_un = atanh( upcorr[`i'] ) - atanh( downcorr[`i'] )
		scalar sigma_Z_un = sqrt(1/( upN[`i'] -3)+1/( downN[`i'] -3))
		qui replace Z_un= mu_Z_un/sigma_Z_un in `i'
		qui replace pvalue_un= 2*normal(-abs(Z_un)) in `i'
		* display "Z statistic = " %8.4g Z _n "P-value = " %8.4g pvalue
	}

keep mypan upcorr upN downcorr downN rawcorr rawN condcorr condN Z pvalue Z_un pvalue_un
keep if !missing(Z)
order mypan upcorr upN downcorr downN rawcorr rawN condcorr condN Z pvalue Z_un pvalue_un

	
* Repeat the test for the diversified factor  

use "anomalies.dta", clear

quietly drop if (ano=="umd" | ano=="glumd" )

qui egen Pan=group(global anomaly)
qui xtset Pan time
qui set more off
qui by Pan,sort:gen n=_n
qui xtset
tssmooth ma MA=ret ,window(12)
drop if n<13
gen CondRet=MA*ret 
rename ret AnomalyRet
merge m:1 year month using "FactorUMD.dta", nogenerate keep(3)

by year month,sort:egen Passive=mean(AnomalyRet)
by year month,sort:egen tempup=mean(AnomalyRet) if MA>0
by year month,sort:egen PassiveUp=mean(tempup) 
by year month,sort:egen tempdown=mean(AnomalyRet) if MA<0
by year month,sort:egen PassiveDown=mean(tempdown) 
drop tempup tempdown
by year month,sort:egen TSMOM=mean(CondRet)

collapse (mean) umd Passive* TSMOM,by(year month)

/* Fisher Test */
corr umd PassiveDown
local i=r(rho)
corr umd PassiveUp
local j=r(rho)
scalar mu_Z = atanh( `j') - atanh( `i')
scalar sigma_Z = sqrt(1/( 666-3)+1/(646-3))
scalar Z = mu_Z/sigma_Z 
scalar pvalue = 2*normal(-abs(Z)) 
display "Z statistic = " %8.4g Z _n "P-value = " %8.4g pvalue


	
*===============================================================================	
* Table IX. Factor Momentum in Momentum-Neutral Factors	
*===============================================================================	

use "oos_tsmom_scs.dta" ,clear

keep yyyymm tsmom1

rename tsmom1 pcmom1

* merge with factor momentum in PCs of momentum-neutral factors
merge 1:1 yyyymm using "oos_tsmom_mom.dta" ,nogenerate
rename tsmom1 pcmom1_neutral
drop tsmom0 tsmom2-tsmom5

foreach var of varlist pcmom1 pcmom1_neutral {
replace `var'=`var'*100
}

qui merge m:1 yyyymm using "FactorFF5.dta", keep(1 3) keepusing(mktrf smb hml rmw cma) nogenerate 

drop if yyyymm<=197306

estimates clear

reg pcmom1 mktrf smb hml rmw cma
estimates store reg5
reg pcmom1 mktrf smb hml rmw cma pcmom1_neutral
estimates store reg6

reg pcmom1_neutral mktrf smb hml rmw cma
estimates store reg7
reg pcmom1_neutral mktrf smb hml rmw cma pcmom1
estimates store reg8

esttab reg*, b(4) t(3) scalars(N r2_a) sfmt(0 5) noparentheses nogaps nostar
	

*-----------------------------------------------
* Table AIII. Conditional Covariances with the Momentum Factor: Decomposition 
*-----------------------------------------------

use "LongShort.dta",clear

egen s_size= rowmean(s_hml s_inv s_op)
egen l_size= rowmean(l_hml l_inv l_op)

reshape long s_ l_,i(year month) j(anomaly) string

replace anomaly="cma" if anomaly=="inv"
replace anomaly="rmw" if anomaly=="op"
replace anomaly="smb" if anomaly=="size"
replace anomaly="strev" if anomaly=="shrev"

egen pan=group(anomaly)

rename l_ retlong
rename s_ retshort

* compute factor return using the existing legs
gen ret=retlong-retshort

gen time=ym(year,month)

qui xtset pan time
qui by pan,sort:gen n=_n
tssmooth ma pastfactorreturn=ret ,window(12)
drop if n<13

merge m:1 year month using "FactorUMD.dta", nogenerate keep(3)

* conditional correlation of factors with UMD 
corr umd ret if pastfactorreturn>0,cov
corr up retlong if pastfactorreturn>0,cov
corr down retshort if pastfactorreturn>0,cov
corr up retshort if pastfactorreturn>0,cov
corr down retlong if pastfactorreturn>0,cov

corr umd ret if pastfactorreturn<0,cov

corr up retlong if pastfactorreturn<0,cov
corr down retshort if pastfactorreturn<0,cov
corr up retshort if pastfactorreturn<0,cov
corr down retlong if pastfactorreturn<0,cov

gen mypan=.

gen positivecorr_umd_ret=.
gen positivecorr_up_retlong=.
gen positivecorr_down_retshort=.
gen positivecorr_up_retshort=.
gen positivecorr_down_retlong=.

gen negativecorr_umd_ret=.
gen negativecorr_up_retlong=.
gen negativecorr_down_retshort=.
gen negativecorr_up_retshort=.
gen negativecorr_down_retlong=.

forvalues i=1/13 {

	* panel 13 is pooled
	if `i'==13 {
		qui replace mypan`char'=`i' in `i'
		
		qui corr  umd ret if pastfactorreturn>0,cov
		qui replace positivecorr_umd_ret=r(cov_12) in `i'
		
		corr up retlong if pastfactorreturn>0,cov
		qui replace positivecorr_up_retlong=r(cov_12) in `i'
		
		corr down retshort if pastfactorreturn>0,cov
		qui replace positivecorr_down_retshort=r(cov_12) in `i'				
		
		corr up retshort if pastfactorreturn>0,cov
		qui replace positivecorr_up_retshort=r(cov_12) in `i'
		
		corr down retlong if pastfactorreturn>0,cov
		qui replace positivecorr_down_retlong=r(cov_12) in `i'
		
		
		qui corr  umd ret if pastfactorreturn<0,cov
		qui replace negativecorr_umd_ret=r(cov_12) in `i'
		
		corr up retlong if pastfactorreturn<0,cov
		qui replace negativecorr_up_retlong=r(cov_12) in `i'
		
		corr down retshort if pastfactorreturn<0,cov
		qui replace negativecorr_down_retshort=r(cov_12) in `i'				
		
		corr up retshort if pastfactorreturn<0,cov
		qui replace negativecorr_up_retshort=r(cov_12) in `i'
		
		corr down retlong if pastfactorreturn<0,cov
		qui replace negativecorr_down_retlong=r(cov_12) in `i'
	}
	else {
		qui replace mypan=`i' in `i'

		qui corr  umd ret if pan==`i' &   pastfactorreturn>0,cov
		qui replace positivecorr_umd_ret=r(cov_12) in `i'
		
		corr up retlong if pan==`i' &  pastfactorreturn>0,cov
		qui replace positivecorr_up_retlong=r(cov_12) in `i'
		
		corr down retshort if pan==`i' &  pastfactorreturn>0,cov
		qui replace positivecorr_down_retshort=r(cov_12) in `i'				
		
		corr up retshort if pan==`i' &  pastfactorreturn>0,cov
		qui replace positivecorr_up_retshort=r(cov_12) in `i'
		
		corr down retlong if pan==`i' &  pastfactorreturn>0,cov
		qui replace positivecorr_down_retlong=r(cov_12) in `i'
		
		
		qui corr  umd ret if pan==`i' &  pastfactorreturn<0,cov
		qui replace negativecorr_umd_ret=r(cov_12) in `i'
		
		corr up retlong if pan==`i' &  pastfactorreturn<0,cov
		qui replace negativecorr_up_retlong=r(cov_12) in `i'
		
		corr down retshort if pan==`i' &  pastfactorreturn<0,cov
		qui replace negativecorr_down_retshort=r(cov_12) in `i'				
		
		corr up retshort if pan==`i' &  pastfactorreturn<0,cov
		qui replace negativecorr_up_retshort=r(cov_12) in `i'
		
		corr down retlong if pan==`i' &  pastfactorreturn<0,cov
		qui replace negativecorr_down_retlong=r(cov_12) in `i'	
	}
}

keep mypan-negativecorr_down_retlong
keep if !missing(negativecorr_down_retlong)
	

*===============================================================================
* Internet appendix: Table IA.I. and Figure IA.1.  factor momentum strategies in off-the-shelf factors
*===============================================================================
use "anomalies.dta", clear

* exclude momentum
drop if (ano=="umd" | ano=="glumd" )

egen Pan=group(global anomaly)
xtset Pan time
by Pan,sort:gen n=_n

set more off

tssmooth ma MA=ret ,window(12)
drop if n<13

* Equal weighted portfolio of all factors
by time,sort:egen EqStrategy=mean(ret)

* TSMOM 
by time,sort:egen TSMomW=mean(ret) if sign(MA)==1
by time,sort:egen TSMomL=mean(ret) if sign(MA)==-1
gen tsmomret=sign(MA)*ret
by time,sort:egen TSMom=mean(tsmomret)

*XSMOM 
* breakpoint for cross-sectional strategy is median 
by time,sort:egen MedianMA=median(MA)
gen signmedian=sign(MA-MedianMA)
replace signmedian=1 if signmedian==0
by time,sort:egen XSMomW=mean(ret) if sign(MA-MedianMA)==1
by time,sort:egen XSMomL=mean(ret) if sign(MA-MedianMA)==-1
gen xsmomret=sign(MA-MedianMA)*ret
by time,sort:egen XSMom=mean(xsmomret)
		
collapse (mean) Eq* TS* XS* ,by(year month)

*Table IA.I.
preserve		
foreach var of varlist Eq* TS* XS* {
	rename `var' ret`var'
}
reshape long ret,i(year month) j(strategy) string
collapse (mean) mean=ret (sd) sd=ret (semean) se=ret (count) n=ret,by(strategy)
gen tvalue=mean/se
replace mean=mean*12
replace sd=sqrt(12)*sd
gen sharpe=mean/sd
list strategy mean sd se tvalue sharpe,sep(0)
* copy to excel
restore	

* Figure IA.1		
* scale strategies such that vol of all strategies is equal to that of the passive strategy 
foreach var of varlist TS* XS* {
	qui sum EqStrategy
	local basevol=r(sd)
	qui sum `var'
	qui replace `var'=`var'*`basevol'/r(sd)
}

foreach	var of varlist EqStrategy TSMomW TSMomL TSMom XSMomW XSMomL XSMom {
	generate sumlogr = sum(ln(1 + `var'/100))
	generate CumRet`var'= exp(sumlogr)  
	* before revision: generate Cumret`var'= exp(sumlogr)
	drop sumlogr
}
	
* OUTPUT for plotting in Matlab
keep year month CumRet*
order year month CumRetEq CumRetTSMomW CumRetTSMomL CumRetXSMomW CumRetXSMomL
* before revision: order year month CumRetEq CumRetTSMomP CumRetTSMomN CumRetXSMomW CumRetXSMomL
sort year month
outsheet using "Temp/timeseries.csv", replace comma


*===============================================================================
*Internet appendix Table IA.II. LM decomposition of TS and XS factor momentum profits
*===============================================================================
 
 /* The TS decomposition */

 
 * Revision: added this paragraph to deal with file not exist problem
 * Check if the file Tests/tmp.dta exists and create if not
capture confirm file Tests/tmp.dta
if _rc {
    clear
    set obs 1  // create a dataset with 1 observation (placeholder)
    gen WAutoCov = .
    gen WMeanSquared = .
    gen bs = .
    save Tests/tmp.dta, replace
}


forvalues bs=0/100 {

	use "anomalies.dta", clear
			
	quietly drop if (ano=="umd" | ano=="glumd" )
	egen Pan=group(anomaly)
	qui xtset Pan time
	set more off
	by Pan,sort:gen n=_n
	qui xtset
	
	qui tssmooth ma MA=ret,window(12)
	qui replace MA=. if n<13
	qui drop if n<13
	drop n

	if `bs'>0 {
		bsample, cluster(time)
	}
	
	/* 1. Autocovariances */
	gen AutoCov=.
	gen N=.
	qui sum Pan
	local x=r(max)
	forvalues j=1/`x'{
		quietly corr ret MA if Pan==`j',covariance
		quietly replace AutoCov`i'=r(cov_12) in `j'
		quietly replace N`i'=r(N) in `j'
	}

	/* 2. Mean squared */
	gen MeanSquared=.
	qui sum Pan
	local x=r(max)
	forvalues j=1/`x'{
		quietly sum ret if Pan==`j'
		quietly replace MeanSquared`i'=r(mean)^2 in `j'
	}

	* weight-adjust contribution by time
	egen WAutoCov= wtmean(AutoCov), weight(N)
	egen WMeanSquared= wtmean(MeanSquared), weight(N) 
	* revision: ssc install _gwtmean
	* this uses David Kantor's _GWTMEAN: Stata module containing extensions to generate to implement weighted mean
	collapse (mean) WAutoCov WMeanSquared 
	gen bs = `bs'
	if `bs'>0 {
		qui append using "Tests/tmp.dta"
	}
	qui save "Temp/tmp.dta", replace 
	
}	
	
	
	
 * The XS decomposition 
 
 * Revision: added this paragraph to deal with file not exist problem
 * Check if the file Temp/tmpcsmean.dta exists and create if not
capture confirm file Temp/tmpcsmean.dta
if _rc {
    clear
    set obs 1  // create a dataset with 1 observation (placeholder)
    gen VarMean = .
    gen bs = .
    save Temp/tmpcsmean.dta, replace
}


 * Auto is the same, compute crosscov and crosssectional variance of premia
 forvalues bs=0/100 {
 
	use "anomalies.dta", clear
			
	quietly drop if (ano=="umd" | ano=="glumd" )
	
	qui egen Pan=group(anomaly)
	qui xtset Pan time
	qui set more off
	qui by Pan,sort:gen n=_n
	qui xtset
	qui tssmooth ma MA=ret,window(12)
	qui replace MA=. if n<13
	qui drop if n<13
	qui drop n

	if `bs'>0 {
		bsample, cluster(time)
	}
	
	/* 1. crosssectional variance of means */
	gen Mean=.
	gen VarMean=.
	forvalues i=1/20{ 
		qui sum ret if Pan==`i'
		qui replace Mean=r(mean) in `i'
	}
	qui sum Mean,d
	replace VarMean=r(Var)

	* COLLAPSE to get the mean
	collapse (mean) VarMean
	gen bs = `bs'
	if `bs'>0 {
		qui append using "Temp/tmpcsmean.dta"
	}
	qui save "Tests/tmpcsmean.dta", replace 
	
}

/* crosscovariances, compute the covariance matrix, remove the diagonal terms, and average the elements */
 forvalues bs=0/100 {
 
	use "anomalies.dta", clear
			
	quietly drop if ( ano=="umd" | ano=="glumd" )
	
	qui egen Pan=group(anomaly)
	qui xtset Pan time
	qui set more off
	qui by Pan,sort:gen n=_n
	qui xtset
	qui tssmooth ma MA=ret,window(12)
	qui replace MA=. if n<13
	qui drop if n<13
	qui drop n
	qui by Pan,sort:gen N=_N

	/* total number of observations are 11,375 */
	
	/* auto and crosscovariances */
	qui keep time ret Pan MA N
	qui reshape wide MA ret N,i(time) j(Pan )
	qui set more off

	if `bs'>0 {
		bsample, cluster(time)
	}

	forvalues i=1/20{
		local obs=1
		qui gen CrossCov`i'=.
		forvalues j=1/20{
			qui corr ret`i' MA`j' ,covariance
			qui replace CrossCov`i'=r(cov_12) in `obs' 
			qui local obs=`obs'+1
		}
	}

	/* The above is the entire matrix. Remove the diagonal to get the cross covs*/

	forvalues i=1/20{
		local obs=1
		forvalues j=1/20{
			qui replace CrossCov`i'=. if `i'==`j' in `obs'
			local obs=`obs'+1
		}
	}
		
	/* scale the crosscovariances by the number of months a factor shows up in the sample,
	this helps to make the summation of components closer to the real premium of the XS strategy
	*/

	gen crosscov=.
	forvalues i=1/20 {
		qui sum CrossCov`i'
		qui replace crosscov=r(mean) in `i'
		qui sum N`i'
		qui replace crosscov=crosscov*r(mean)/11375 in  `i'
	}
		
	collapse (sum) crosscov 
	
	gen bs = `bs'
	if `bs'>0 {
		qui append using "Temp/tmpcscross.dta"
	}
	qui save "Temp/tmpcscross.dta", replace 
} 
 
 
use "Temp/tmp.dta", clear
merge 1:1 bs using "Temp/tmpcscross.dta", nogenerate
merge 1:1 bs using "Temp/tmpcsmean.dta", nogenerate

* Total XS strategy (= autovariance + crosscovariance + cross-sectional variance)

use "anomalies.dta", clear
quietly drop if (ano=="umd" | ano=="glumd" )
	
qui egen Pan=group(anomaly)
qui xtset Pan time
qui set more off
qui by Pan,sort:gen n=_n
qui xtset

tssmooth ma MA=ret,window(12)
replace MA=. if n<13
drop if n<13
drop n

by time,sort:egen AveMA=mean(MA)
by time,sort:gen  NoAllStrategies=_N
gen Weight=1/NoAllStrategies*(MA-AveMA)
by time,sort:egen XSRet=total(Weight*ret)
collapse (mean) XSRet,by(year month)

reg XSRet


*===============================================================================
* Internet appendix Table IA.III. 
*===============================================================================

* create blank files for the TS and XS strategy results
clear
gen HoldingPeriod = .
save "Temp/FormHoldTS.dta", replace
save "Temp/FormHoldXS.dta", replace

* {h,k} formation/holding for Time Series Strategy 
forvalues i=1/24 {

	use "anomalies.dta", clear
	qui drop if (ano=="umd" | ano=="glumd" )
	qui egen Pan=group(global anomaly)
	qui xtset Pan time
	qui by Pan,sort:gen n=_n
	qui set more off

	qui tssmooth ma MA`i' =  ret, window(`i')
	qui replace MA`i'=. if n<`i'
	qui gen flag=1 if f.MA`i'>0 & f.MA`i'!=.
	qui replace flag=0 if (f.MA`i'<0 | f.MA`i'==0) & f.MA`i'!=.
	qui drop if f.MA`i'==.

	* Find the `Active' portfolios 
	forvalues j=1/24{
		qui xtset
		qui gen ActiveRet`j'=ret if l`j'.flag==1
		qui replace ActiveRet`j'=-1*ret if l`j'.flag==0
		qui by time,sort:egen AveActive`j'=mean(ActiveRet`j')
	}

	qui drop if time==time[_n-1]

	gen alpha`i'=.
	gen tstat`i'=.
	quietly tsset time
	forvalues x=1/24{
		quietly reg AveActive`x' 
		qui replace alpha`i'=_b[_cons] in `x'
		qui replace tstat`i'=_b[_cons]/_se[_cons] in `x'
	}
	qui keep alpha* tstat*
	qui gen HoldingPeriod=_n
	qui merge 1:1 HoldingPeriod using "Temp/FormHoldTS.dta"
	qui drop _m
	qui drop if HoldingPeriod>24
	save "Temp/FormHoldTS.dta",replace

}

order HoldingPeriod alpha* tstat*, sequential
format alpha* tstat* %12.2f

**************************************************************************
* {h,k} formation/holding for the cross sectional strategy 

forvalues i=1/24 {
	
	use "anomalies.dta", clear
	quietly drop if (ano=="umd" | ano=="glumd" )
	quietly egen Pan=group(global anomaly)
	quietly xtset Pan time
	quietly by Pan,sort:gen n=_n
	quietly set more off
	
	* compute average return over the past `i' months
	qui tssmooth ma MA`i' =  ret, window(`i')
	qui replace MA`i'=. if n<=`i'
	
	* compute the median of average returns across all anomalies
	qui by time,sort:egen medMA`i'=median(MA`i')
	qui xtset
	
	* +1/-1 strategy 
	qui gen flag=sign(f.MA`i'-f.medMA`i')
	
	* all anomalies are in either in the long or the short leg
	replace flag=1 if flag==0
	
	* Find the `Active' portfolios 
	forvalues j=1/24 {
		qui xtset
		qui gen ActiveRet`j'=ret*l`j'.flag
		quietly by time,sort:egen AveActive`j'=mean(ActiveRet`j')
	}
	qui collapse (mean) AveActive*,by(time)
		
	qui gen alpha`i'=.
	qui gen tstat`i'=.
	qui tsset time
	qui forvalues x=1/24 {
		qui reg AveActive`x' 
		qui replace alpha`i'=_b[_cons] in `x'
		qui replace tstat`i'=_b[_cons]/_se[_cons] in `x'
	}
	
	qui keep alpha* tstat*
	qui gen HoldingPeriod=_n
	qui merge 1:1 HoldingPeriod using "Temp/FormHoldXS.dta", nogenerate
	qui drop if HoldingPeriod>24
	qui save "Temp/FormHoldXS.dta", replace
	
}

order HoldingPeriod alpha* tstat*, sequential
format alpha* tstat* %12.2f


*===============================================================================
* internet appendix IA.III. Decomposing Equity Momentum Profits
*
* The input file is the monthly CRSP file from WRDS
* The paper uses the file that ended in December 2019
*===============================================================================
* construct 100 momentum sorted portfolios

/*
use "Monthly1926-2019Permco.dta", clear
keep if (shrcd==10 | shrcd==11)
gen time=ym(year,month)
drop if time<42
xtset permco time
bys permco: asrol ret, stat(mean) win(time 11) min(11)
gen mom=l2.mean11_ret
gen me=abs(prc* shrout)
gen l_me=l.me
drop if (missing(mom) | missing(l_me) | missing(ret))
keep permco exchcd l_me mom year month time ret
save "Temp/FirmLevelMom.dta", replace

* 100 portfolios sorted on momentum (pre-ranking and post-ranking returns)
use "Temp/FirmLevelMom.dta", clear
gen nyse_SORTVAR = mom if exchcd==1
forvalues p=1/99 {
	bysort time: egen pct`p' = pctile(nyse_SORTVAR), p(`p')
}

drop nyse_SORTVAR

gen long Percentile = 1 if mom <= pct1 & ~missing(pct1)
forvalues p=2/99 {
	local pp = `p' - 1
	qui replace Percentile = `p' if (mom > pct`pp' & mom <= pct`p' & ~missing(pct`pp') & ~missing(pct`p'))
} 

replace Percentile = 100 if mom > pct99 & ~missing(pct99) & ~missing(mom)
drop if missing(Percentile)
	
egen PastRet= wtmean(mom), by(Percentile year month) weight(l_me)
egen Ret= wtmean(ret), by(Percentile year month) weight(l_me)
collapse (mean) PastRet Ret, by(Percentile year month)
replace PastRet=(PastRet)*100
replace Ret=Ret*100
save "Temp/Decomposition_Returns100.dta", replace

* implement LM
bys year month: egen AvePastRet=mean(PastRet)
bys year month: egen AveRet=mean(Ret)
gen time=ym(year,month)
gen Strategy = (PastRet - AvePastRet) * (Ret)
collapse (mean) Strategy, by(year month)
save "Temp/Decomposition_Strategy100.dta", replace

* compute past factor returns 
use "FactorFF5.dta", clear
merge 1:1 year month using "FactorBAB.dta", keep(match) nogenerate
merge 1:1 year month using "FactorQMJ.dta", keep(match) nogenerate
gen time=ym(year,month)
tsset time
foreach var of varlist mktrf-cma bab qmj {
	qui tssmooth ma temp`var'=`var',window(11)
	gen MA`var'=l.temp`var'
}

keep year month MA* 
save "Temp/PastFF5Ret.dta",replace
 
* portfolio daily data for estimation of betas 
use "Temp/FirmLevelMom.dta", clear

gen nyse_SORTVAR = mom if exchcd==1
forvalues p=1/99 {
	bysort time: egen pct`p' = pctile(nyse_SORTVAR), p(`p')
}
drop nyse_SORTVAR

gen long Percentile = 1 if mom <= pct1 & ~missing(pct1)
forvalues p=2/99 {
	local pp = `p' - 1
	qui replace Percentile = `p' if (mom > pct`pp' & mom <= pct`p' & ~missing(pct`pp') & ~missing(pct`p'))
} 

replace Percentile = 100 if mom > pct99 & ~missing(pct99) & ~missing(mom)
drop if missing(Percentile)
keep year month permco l_me Percentile
merge 1:m permco year month  using "Daily1926-2019_Permco.dta", nogenerate
egen PortDailyRet= wtmean(ret), by(Percentile year month date) weight(l_me)
gen day=day(date)
collapse (mean) PortDailyRet, by(Percentile year month day)
gen time=ym(year,month)
replace PortDailyRet= PortDailyRet*100
save "Temp/DailyReturnsPercentile.dta", replace

* estimate betas
use "Temp/DailyReturnsPercentile.dta", clear

merge m:1 year month day using "DailyFactorFF5.dta", keep(match) nogenerate
merge m:1 year month day using "DailyFactorBAB.dta", keep(match) nogenerate
merge m:1 year month day using "DailyFactorQMJ.dta", keep(match) nogenerate
gen retrf=PortDailyRet-rf

bys Percentile: asreg retrf mktrf , wind(time 3) 
rename _b_mktrf b_mktrf
drop _*

bys Percentile: asreg retrf mktrf smb hml, wind(time 3)  
rename _b_mktrf b_mktrfFF3
rename _b_smb b_smbFF3
rename _b_hml b_hmlFF3
drop _*

bys Percentile: asreg retrf mktrf smb hml rmw cma, wind(time 3)  
rename _b_mktrf b_mktrfFF5
rename _b_smb b_smbFF5
rename _b_hml b_hmlFF5
rename _b_rmw b_rmwFF5
rename _b_cma b_cmaFF5
drop _*

bys Percentile: asreg retrf mktrf smb hml rmw cma bab qmj, wind(time 3) 
rename _b_mktrf b_mktrfFFBAB
rename _b_smb b_smbFFBAB
rename _b_hml b_hmlFFBAB
rename _b_rmw b_rmwFFBAB
rename _b_cma b_cmaFFBAB
rename _b_bab b_babFFBAB
rename _b_qmj b_qmjFFBAB
drop _*
drop bab
collapse (mean) b* ,by(Percentile year month)
	
save "Temp/MonthlyCoefficientsPercentile.dta", replace

*/


**************************************************************
* Estimate the matrix using 100 momentum sorted portfolios 
* compute standard errors for each component using bootstraps
* The residual component is defined as (premium - sum of other components)



* Revision: added this code block to check if the file Temp/tmp_bs_decomposition.dta exists and create if not

capture confirm file Temp/tmp_bs_decomposition.dta
if _rc {
	clear
	set obs 1  // create a dataset with 1 observation (placeholder)
	gen bs = .
	save Temp/tmp_bs_decomposition.dta, replace
}


forvalues bs=0/1000 {

	qui use "Temp/MonthlyCoefficientsPercentile.dta", clear
	qui merge 1:1 Percentile year month using "Temp/Decomposition_Returns100.dta", keep(match) nogenerate
	qui merge m:1 year month using "Temp/PastFF5Ret.dta", keep(match) nogenerate
	qui merge m:1 year month using "FactorFF5.dta", keep(match) nogenerate
	qui merge m:1 year month using "FactorBAB.dta", keep(match) nogenerate
	qui merge m:1 year month using "FactorQMJ.dta", keep(match) nogenerate
	qui merge m:1 year month using "Temp/Decomposition_Strategy100.dta", keep(match) nogenerate
	qui gen time=ym(year,month)
	
	if `bs'>0 {
		bsample, cluster(time)
	}
		
	************ 7 Factor model *************************************** 

	* covariance matrix of betas
	qui corr b_mktrfFFBAB b_smbFFBAB b_hmlFFBAB b_rmwFFBAB b_cmaFFBAB b_babFFBAB b_qmjFFBAB,cov
	matrix BetaF7=r(C)

	*   factor cross-covariances 
	matrix FacCov = J(7,7,0)
	qui matrix list FacCov, nohalf

	local varcount=1
	qui foreach var of varlist MAmktrf MAsmb MAhml MArmw MAcma MAbab MAqmj {
		qui corr `var' mktrf smb hml rmw cma bab qmj,cov
		qui matrix tempcov=r(C)
		forvalues i=1/7 {
			qui matrix FacCov[`i',`varcount']=tempcov[`i'+1,1]
		}

		local varcount=`varcount'+1
	}
		
	* compute the sum of all products
	matrix sumcomp1comp2= J(7,7,0)
	forvalues j=1/7 {
		forvalues i=1/7 {
			qui matrix sumcomp1comp2[`i',`j']=BetaF7[`i',`j']*FacCov[`i',`j']
		}
	}
		
	qui mat li sumcomp1comp2, nohalf

	local comp1 = trace(sumcomp1comp2)
	qui gen comp1_7F=`comp1'

	qui matsum sumcomp1comp2,all(sumscalar) display
	local sum=sumscalar[1,1]
	qui gen comp2_7F=`sum'-`comp1'

	clear matrix
		
	************ 5 Factor model *************************************** 

	* covariance matrix of betas
	qui corr b_mktrfFF5 b_smbFF5 b_hmlFF5 b_rmwFF5 b_cmaFF5 ,cov
	qui matrix BetaF5=r(C)

	*  factor cross-covariances 
	qui matrix FacCov = J(5,5,0)

	local varcount=1
	qui foreach var of varlist MAmktrf MAsmb MAhml MArmw MAcma {
		qui corr `var' mktrf smb hml rmw cma ,cov
		qui matrix tempcov=r(C)
		forvalues i=1/5 {
			qui matrix FacCov[`i',`varcount']=tempcov[`i'+1,1]
		}
		local varcount=`varcount'+1
	}
		
	* compute the sum of all products
	matrix sumcomp1comp2= J(5,5,0)
	forvalues j=1/5 {
		forvalues i=1/5 {
			qui matrix sumcomp1comp2[`i',`j']=BetaF5[`i',`j']*FacCov[`i',`j']
		}
	}
		
	qui mat li sumcomp1comp2, nohalf

	local comp1 = trace(sumcomp1comp2)
	qui gen comp1_5F=`comp1'

	qui matsum sumcomp1comp2,all(sumscalar) display
	local sum=sumscalar[1,1]
	qui gen comp2_5F=`sum'-`comp1'

	clear matrix
		
	************ 3 Factor model *************************************** 

	* covariance matrix of betas
	qui corr b_mktrfFF3 b_smbFF3 b_hmlFF3,cov
	matrix BetaF3=r(C)

	*  factor cross-covariances 
	matrix FacCov = J(3,3,0)
	qui matrix list FacCov , nohalf

	local varcount=1
	qui foreach var of varlist MAmktrf MAsmb MAhml {
		qui corr `var' mktrf smb hml ,cov
		qui matrix tempcov=r(C)
		forvalues i=1/3 {
			qui matrix FacCov[`i',`varcount']=tempcov[`i'+1,1]
		}
		local varcount=`varcount'+1
	}
		
	* compute the sum of all products
	matrix sumcomp1comp2= J(3,3,0)
	forvalues j=1/3 {
		forvalues i=1/3 {
			qui matrix sumcomp1comp2[`i',`j']=BetaF3[`i',`j']*FacCov[`i',`j']
		}
	}
		
	qui mat li sumcomp1comp2, nohalf

	local comp1 = trace(sumcomp1comp2)
	qui gen comp1_3F=`comp1'

	qui matsum sumcomp1comp2,all(sumscalar) display
	local sum=sumscalar[1,1]
	qui gen comp2_3F=`sum'-`comp1'

	clear matrix
		
	************ CAPM model *************************************** 

	* covariance matrix of betas
	qui corr b_mktrf,cov
	matrix BetaF1=r(C)

	*   factor cross-covariances 
	qui matrix FacCov = J(1,1,0)
	qui matrix list FacCov , nohalf

	local varcount=1
	qui foreach var of varlist MAmktrf {
		qui corr `var' mktrf ,cov
		qui matrix tempcov=r(C)
		forvalues i=1/1 {
			qui matrix FacCov[`i',`varcount']=tempcov[`i'+1,1]
		}
		local varcount=`varcount'+1
	}
		
	* compute the sum of all products
	qui matrix sumcomp1comp2= J(1,1,0)
	forvalues j=1/1 {
		forvalues i=1/1 {
			qui matrix sumcomp1comp2[`i',`j']=BetaF1[`i',`j']*FacCov[`i',`j']
		}
	}
		
	local comp1 = trace(sumcomp1comp2)
	qui gen comp1_1F=`comp1'

	qui matsum sumcomp1comp2,all(sumscalar) display
	local sum=sumscalar[1,1]
	qui gen comp2_1F=`sum'-`comp1'

	clear matrix

	bys Percentile: egen meanret=mean(Ret)

	qui collapse (mean) comp* meanret Strategy,by(Percentile)

	qui collapse (mean) comp* Strategy (sd) meanret
	qui gen comp4=meanret^2
	qui gen sumcomp1comp2_1f=comp1_1F+ comp2_1F
	qui gen sumcomp1comp2_3f=comp1_3F+ comp2_3F
	qui gen sumcomp1comp2_5f=comp1_5F+ comp2_5F
	qui gen sumcomp1comp2_7f=comp1_7F+ comp2_7F

	gen comp3_7F=Strategy-comp1_7F-comp2_7F-comp4
	gen comp3_5F=Strategy-comp1_5F-comp2_5F-comp4
	gen comp3_3F=Strategy-comp1_3F-comp2_3F-comp4
	gen comp3_1F=Strategy-comp1_1F-comp2_1F-comp4

	keep comp* sum* 

	order comp1* comp2* comp3* comp4* 
	

	qui gen bs = `bs'
	if `bs'>0 {
		qui	append using "Temp/tmp_bs_decomposition.dta"
	}
	qui	save "Temp/tmp_bs_decomposition.dta", replace
	
}
