*{ 	//############# PIVUS ##############(

*{	//#############  prepare data for analysis ##############(
local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
use "Endostatin PIVUS 10 year follow-up.dta", clear 
drop _merge

duplicates drop id, force
merge 1:1 id using "PIVUS70 prevalent svikt droppadeStart.dta"
drop _merge
save "Endostatin PIVUS 10 year follow-up merged PIVUS70 .dta", replace

insheet using "pivus_olink_qc_notransform_LODreplace.txt", clear 
merge 1:1 id using "Endostatin PIVUS 10 year follow-up merged PIVUS70 .dta"

histogram allmidatestata 
tab allmi if beskdat70mdy < allmidatestata

tab hjrtsvikt svikt // prevalent/history of heart failure at baseline 
tab svikt // svikt = heart failure during follow-up 

drop if hjrtsvikt ==1 // remove prevalent heart failure 
drop if prot_il_8 ==. // remove individuals without Proseek Multiplex CVD I measurement 
drop prot
drop protein
keep if _merge==3
drop _merge
// drop proteins missing in either cohort 
	drop prot_hsp_27
	drop prot_cstb

g waist=mage 
g male=0
recode male 0=1 if kn==0 
label var male "Male=1"
g tg=triglycerider
g age=exaktlder70
g lipmed =0
replace lipmed=1 if statiner==1
replace lipmed=1 if andralipidmed==1
replace afib=1 if AFib==1
generate gfr=(77.24*cystc^-1.2623)
	// GFR in mL/min can be calculated from serum cystatin C results in mg/L by the formula y = 77.24x-1.2623. The assays were performed at the Department of Clinical Chemistry, University Hospital, Uppsala, which is accredited according to 17025 by Swedac.
	// cystc: SERUM CYSTATIN C (mg/L)

rename insulin insulintreatment
rename rkarenu smoker
egen stdmage=std(waist) 
egen stdntprobnp=std(lnprobnp) 
g prevAF=afib 

gen non_isc_hf = allsvikt 
recode non_isc_hf 1=0 if validmi70 == 1 
recode non_isc_hf 1=0 if allmidatestata <= sviktdatestata 
gen isc_hf= allsvikt 
recode isc_hf 1=0 if non_isc_hf == 1 
g fastingglucose=(blodsocker*1.11) 
gen y77 = real(substr(deathdate,1,4)) 
gen m77 = real(substr(deathdate,6,2))
gen d77 = real(substr(deathdate,9,2))
gen deathdatenew = mdy(m77, d77, y77) 
format deathdatenew %d
label variable deathdatenew "date of death in Stata format"
drop m77 d77 y77 
sum id 

g allsviktdatestataoriginal = allsviktdatestata 

//removing missing values and making exclusion table 
// removed waist from exclusion due to many missing in ULSAM 

foreach var in lipmed medicinmothgtbt oralaantidiabetika insulintreatment smoker { // These have been verified by Johan to be miscoded 
	replace `var' = 0 if `var' ==. 
	}

file open table2 using "MissingDataPIVUS.txt", write replace
foreach var in bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp medicinmothgtbt smoker SV2RV534 gfr ntprobnp {
	sum id if `var' ==. 
	file write table2 "`var'" _tab (r(N)) _n
	drop if `var' ==. 
	}

file close table2  

// Used when updating flowchart then commented out:
	// insheet using "MissingDataPIVUS.txt", clear 
	// insheet using "MissingDataULSAM.txt", clear // For updating flowchart with confounder exclusion 
	// export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("4pExclusions") 

drop if allsviktdatestata <= beskdat70mdy 
// removed id 654 929 974 (dropped because event before start of study)
save "PIVUS70 prevalent svikt droppade Before Stsplit.dta", replace 

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

*{	//##############  make time dependent covariates ##############(
 
*{	//##############  stsplit on myocardial infarction ##############(

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)
recode allmi 1=0 if allmidatestata >= allsviktdatestata // All MI after heart failure coded as 0
recode allmi 0=1 if validmi70==1 

gen amistart=allmidatestata if allmi==1 
format amistart %d
replace amistart=beskdat70mdy if validmi70==1 // validmi70 are all cases of MI that occured on or before start time (renamed to prevMI
replace amistart=beskdat70mdy if allmidatestata<beskdat70mdy // 4 cases where allmidatestata indicate that individuals had MI before start of follow-up 
tab validmi70 allmi // 49 individuals had MI at start of study and 94 (109 including those having MI after HF) at end of study 

sort id 
g start1=beskdat70mdy
format start1 %d
g stop1=amistart if allmi==1 
format stop1 %d
g start2=amistart+1 if allmi==1 
format start2 %d 
g stop2=allsviktdatestata 
format stop2 %d

g mifreetime = (amistart-beskdat70mdy) // shows how many days each individual was free of MI during the study 
replace mifreetime = (allsviktdatestata-beskdat70mdy) if allmi==0 & mifreetime==.

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

replace amistart = 100000 if amistart == . 
stsplit HaveMI, at(0) after(amistart) 
replace HaveMI = HaveMI + 1

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)


*}
	//############End stsplit on myocardial infarction ##############)

*{	//##############  stsplit on atrial fibrillation ##############

g afibstart=afibdatestata if afib==1 
format afibstart %d
replace afibstart=beskdat70mdy if afibdatestata<beskdat70mdy // cases where  individuals had A-fib before start of follow-up  (All afib after heart failure coded as 0)

sum afibstart 
replace afibstart = 100000 if afibstart == . 
stsplit HaveAfib, at(0) after(afibstart) 
replace HaveAfib = HaveAfib + 1

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

duplicates tag id, g(duplicates)
list id _t0 _t _d _st HaveMI HaveAfib amistart afibstart if duplicates == 1 

save "PIVUS70 prevalent svikt droppade.dta", replace  

*}
	//############End stsplit on atrial fibrillation ##############)
 	
*} 
	//############End make time dependent covariates ##############)

*}
	//###########End prepare data for analysis ##############)

*{	//#############  Primary Analysis (

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
use "PIVUS70 prevalent svikt droppade.dta", clear 
stset allsviktdatestata, id(id) fail(allsvikt==1) enter(time beskdat70mdy) origin(time birthdate) scale(365.25) 

file open table2 using "Discovery_PIVUS.txt", write replace  
file write table2 "Variable" ///
		_tab "Lower 95% CI" ///
		_tab "Age Adjusted HR" ///
		_tab "upper 95% CI" ///
		_tab "P-value1"  ///
		_tab "betaP-value1" ///
		_tab "seP-value1" ///
		_tab "N_tot" ///
		_tab "N_events" ///		
	_tab "Variable" ///
		_tab "Multiple Adjusted Lower 95% CI" ///  
		_tab "Multiple Adjusted HR" ///
		_tab "Multiple Adjusted upper 95% CI" ///
		_tab "P-value2" ///
		_tab "beta" ///
		_tab "se" ///
		_tab "Multiple Adjusted N_tot" ///
		_tab "Multiple Adjusted N_events" _n

	foreach var of varlist prot*  { // The following is a Cox-regression for all proteins 
		quietly:stcox `var' kn  
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///			
			_tab
		quietly:stcox `var' HaveMI HaveAfib manuelltsbp kn medicinmothgtbt hdl ldl lipmed smoker diabetesmellitus validmi70 bmi SV2RV534 gfr afib 
		// quietly:stcox `var' manuelltsbp kn medicinmothgtbt hdl ldl lipmed smoker diabetesmellitus validmi70 bmi SV2RV534 gfr afib ntprobnp
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) _n
		}

file close table2 
insheet using "Discovery_PIVUS.txt", clear
g cohort="PIVUS"
multproc, pvalue(pvalue1) method(simes) nhcred(signpvalue1)
multproc, pvalue(pvalue2) method(simes) nhcred(signpvalue2)
sort pvalue1 

sort variable 
save "Regression table_proteomics vs svikt_PIVUS_HR.dta", replace
// preserve 
g pvalue1PIVUS = pvalue1
g signpvalue1PIVUS = signpvalue1
keep variable pvalue1PIVUS signpvalue1 
save "Significant variables PIVUS.dta", replace
// restore 

*}
	//###########End Primary Analysis )  

*} 
	//######### End PIVUS ##############)

*{	//#############  ULSAM ##############(

*{ 	//#############  prepare data for analysis ###########( 

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'

use "ULSAM data proteomics for individual metaanalysis_Afib_korr_2.dta", clear 
rename lpnr id
keep afib afibdatestata id
save "ULSAM data proteomics for individual metaanalysis_Afib_korr_2 keeping afib afibdatestata id.dta", replace 

use "ULSAM77 STATA12", clear 
drop _merge
rename pat id 
save "ULSAM77 STATA12 Temp", replace 
insheet using "ulsam_olink_qc_notransform_LODreplace.txt", clear 

merge 1:1 id using "ULSAM77 STATA12 Temp.dta"
sort id
keep if _merge==3
drop _merge
merge m:m id using "ULSAM data proteomics for individual metaanalysis_Afib_korr_2 keeping afib afibdatestata id.dta" 

drop if _merge==2 // _merge==1 & 3 should be kept (has been confirmed) 
drop _merge
drop if age77==. // remove individuals not participating in ULSAM77 

gen y77 = real(substr(dat77,1,4)) 
gen m77 = real(substr(dat77,5,2))
gen d77 = real(substr(dat77,7,2))
gen sdat77 = mdy(m77, d77, y77)
format sdat77 %d
drop m77 d77 y77 
label variable sdat77 "statadatum baseline77"
gen y70 = real(substr(sviktdat,1,4)) 
gen m70 = real(substr(sviktdat,5,2))
gen d70 = real(substr(sviktdat,7,2))
gen sviktdat1 = mdy(m70, d70, y70)
format sviktdat1 %d
drop m70 d70 y70 
label variable sviktdat1 "statadatum chf"

drop if sdat77>sviktdat1 // remove prevalent heart failure 

gen y77 = real(substr(p064,1,4)) 
gen m77 = real(substr(p064,5,2))
gen d77 = real(substr(p064,7,2))
gen amidate = mdy(m77, d77, y77)
format amidate %d
drop m77 d77 y77 
label variable amidate "statadatum amidate"
gen y70 = real(substr(p064,1,4)) 
gen m70 = real(substr(p064,5,2))
gen d70 = real(substr(p064,7,2))
gen samidat1 = mdy(m70, d70, y70)
format samidat1 %d
drop m70 d70 y70 
label variable samidat1 "statadatum ami"
gen y77 = real(substr(bdate,1,4)) 
gen m77 = real(substr(bdate,5,2))
gen d77 = real(substr(bdate,7,2))
gen birthdate = mdy(m77, d77, y77)
format birthdate %d
drop m77 d77 y77 
label variable birthdate "statadatum birthdat"
gen y77 = real(substr(dat77,1,4)) 
gen m77 = real(substr(dat77,5,2))
gen d77 = real(substr(dat77,7,2))

format sdat77 %d
drop m77 d77 y77 
label variable sdat77 "statadatum baseline77"
gen y70 = real(substr(p013,1,4)) 
gen m70 = real(substr(p013,5,2))
gen d70 = real(substr(p013,7,2))

gen smortdat1 = mdy(m70, d70, y70)
format smortdat1 %d
drop m70 d70 y70 
label variable smortdat1 "statadatum mortalitet"


gen y70 = real(substr(p082,1,4)) 
gen m70 = real(substr(p082,5,2))
gen d70 = real(substr(p082,7,2))
gen sstrodat1 = mdy(m70, d70, y70)
format sstrodat1 %d
drop m70 d70 y70 
label variable sstrodat1 "statadatum stroke"
gen y70 = real(substr(sviktdat,1,4)) 
gen m70 = real(substr(sviktdat,5,2))
gen d70 = real(substr(sviktdat,7,2))

format sviktdat1 %d
drop m70 d70 y70 
label variable sviktdat1 "statadatum chf"

gen allCVD=0
recode allCVD 0=1 if p063==1 
recode allCVD 0=1 if p063==2
recode allCVD 0=1 if p081==1
recode allCVD 0=1 if p081==2
recode allCVD 0=1 if svikt==1
tab allCVD

//Make variable diabetes (1/0)
gen byte diabetes = .
replace diabetes = 0 if pglukos <7.0 & pglukos !=. 
replace diabetes = 1 if pglukos !=. & pglukos >7|pglukos ==7
replace diabetes = 1 if diabmed ==1

gen amiprev = p063
recode amiprev 1=0 if sdat77<samidat1
recode amiprev 2=0
gen strokeprev = p081
recode strokeprev 1=0 if sdat77<sstrodat1 
gen chfprev = svikt
recode chfprev 1=0 if sdat77<sviktdat1 
drop if chfprev ==1  

//generate start time analysis for each patient
gen allsviktdatestata = mdy(12, 31, 2008)
format allsviktdatestata %d
label variable allsviktdatestata "Exit date svikt"

recode p001 1=0 if smortdat1>allsviktdatestata
recode p001 2=0 if smortdat>allsviktdatestata

replace smortdat1 = allsviktdatestata if smortdat1>allsviktdatestata
replace allsviktdatestata = smortdat1 if smortdat1!=.
replace allsviktdatestata = sviktdat1 if sviktdat1!=.
duplicates list id 
duplicates drop 

gen lvh77=0
recode lvh77 0=1 if  v724==1 
recode lvh77 0=1 if  v725x==1
recode lvh77 0=1 if  v725y==1
tab lvh77

gen hyptreat =0
replace hyptreat = 1 if alpha ==1|diur ==1|beta ==1|ca ==1|ace ==1
rename a16 smoker
rename v324 ldl 
rename v302 hdl
rename svikt allsvikt
rename v013 manuelltsbp
rename x1 kolesterol
rename diabetes diabetesmellitus
rename v020 mage
rename v290 bmi
rename sdat77 beskdat70mdy
gen kn=1
rename v971 triglycerider
g manuelltdbp = v014
rename hyptreat medicinmothgtbt
rename amiprev validmi70
rename lvh77 SV2RV534
recode smoker .=2
g tg=triglycerider
rename v972 chol
rename v974 WHR 
rename prot_nt_ ntprobnp 
//gen lipmed =0
recode lipmed 0=1 if statin==1
recode lipmed 0=1 if v105==1
// g fastingglucose = blodsocker 
// g glukos=exp(lnglukos) 
// tab glukos blodsocker 
g fastingglucose = v319 
rename v406 insulintreatment  
rename v407 oralaantidiabetika 
rename v730 prevAF
// g afib=v730
recode afib .=0
recode prevAF .=0
recode prevAF 9=0
// g manuelltdbp = v014
// g age=age77
g waist=mage
generate gfr=(77.24*v800^-1.2623)
	// GFR in mL/min can be calculated from serum cystatin C results in mg/L by the formula y = 77.24x-1.2623. The assays were performed at the Department of Clinical Chemistry, University Hospital, Uppsala, which is accredited according to 17025 by Swedac.
	// V800	77Y: SERUM CYSTATIN C (mg/L)
	// V809	77Y: URINE CYSTATIN C (mg/L)

	egen stdntprobnp=std(ntprobnp)

// drop if missing in either cohort 
drop prot_beta_ngf
drop prot_sirt2
drop prot_nemo
drop ptx3
drop prot_mmp_7

rename age77 exaktlder70 
gen non_isc_hf = allsvikt 
recode non_isc_hf 1=0 if validmi70 == 1 
recode non_isc_hf 1=0 if samidat1 <= sviktdat1
gen isc_hf = allsvikt 
recode isc_hf 1=0 if non_isc_hf == 1

//removing missing values and making exclusion table 
	file open table2 using "MissingDataULSAM.txt", write replace
		foreach var in lipmed medicinmothgtbt oralaantidiabetika insulintreatment insulin diabetesmellitus { // not dropping those that have been verifiedby Johan to be miscoded 
			replace `var' = 0 if `var' ==. 
			}

		foreach var in kn tg hdl ldl smoker validmi70 bmi gfr manuelltsbp manuelltdbp SV2RV534 lipmed fastingglucose insulintreatment oralaantidiabetika medicinmothgtbt ntprobnp {
		sum id if `var' ==. 
		file write table2 "`var'" _tab (r(N)) _n
		drop if `var' ==. 
	// make std(var) 
		// egen std`var'=std(`var') 
		}
	file close table2

save "ULSAM77 data proteomics Proseek after exclusions.dta", replace 

// Used when updating flowchart then commented out 
	// insheet using "MissingDataULSAM.txt", clear
	// export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("4uExclusions") 

drop if allsviktdatestata <= beskdat70mdy 
g allmi = p056 

histogram allmidatestata 
tab allmi if beskdat70mdy < allmidatestata 

save "ULSAM77 data proteomics Proseek Before Stsplit.dta", replace 

*{	//##############  make time dependent covariates ##############(

*{	//##############  stsplit on myocardial infarction ##############(

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)
recode allmi 1=0 if amidate >= allsviktdatestata // To have all MI after heart failure coded as 0. 
recode allmi 0=1 if validmi70==1 

tab validmi70 allmi // 49 individuals had MI at start of study and 94 (including tohse that had MI after HF it is 109) at end of study 
g amistart=amidate if allmi==1 


replace amistart = 100000 if amistart == . 
stsplit HaveMI, at(0) after(amistart) 
replace HaveMI = HaveMI + 1
format amistart %d 

duplicates tag id, g(duplicates)
list id _t0 _t _d _st amistart HaveMI if duplicates == 1

replace allsvikt =0 if allsvikt == .  

stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

stcox HaveMI exaktlder70 manuelltsbp kn medicinmothgtbt hdl ldl lipmed smoker diabetesmellitus validmi70 bmi SV2RV534 gfr afib ntprobnp 
*}
	//############End stsplit on myocardial infarction ##############)

*{	//##############  stsplit on atrial fibrillation ##############(
// v730 are cases of atrial fibrillation at start of study. These individuals have missing afibdatestata
// afib are cases of atrial fibrillation not at start of study (unclear start and stop) 
g allafib=0 
recode allafib 0=1 if afib==1 | prevAF ==1 // allafib includes all atrial fibrillation cases 

g afibduringstudy=allafib // variable for events during study 
recode afibduringstudy 1=0 if afibdatestata>=allsviktdatestata // cases where  individuals had A-fib before start of follow-up (All afib after heart failure coded as 0)

gen afibstart=afibdatestata if allafib==1 
format afibstart %d

replace afibstart = 100000 if afibstart == . | afibdatestata>=allsviktdatestata
replace afibstart=beskdat70mdy if prevAF==1  // cases where  individuals had A-fib before start of follow-up 
replace afibstart=beskdat70mdy if afibdatestata<=beskdat70mdy // x cases where  individuals had A-fib before start of follow-up 
sort allafib prevAF afib afibstart 

stsplit HaveAfib, at(0) after(afibstart) 
replace HaveAfib = HaveAfib + 1

tab allafib HaveAfib
list id _t0 _t _d _st amistart HaveMI HaveAfib if duplicates == 1 
 
stset allsviktdatestata, fail(allsvikt == 1) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

save "ULSAM77 data proteomics Proseek.dta", replace 

*}
	//############End stsplit on atrial fibrillation ##############)
 	
 *} 
	//############End make time dependent covariates ##############)

*{	//##############  Removing all cases of atherosclerosis ##############(
// use "ULSAM77 data proteomics Proseek.dta", clear 
drop if p054 > 0 
drop if p072 > 0 
save "ULSAM77 data proteomics Proseek without atherosclerosis.dta", replace 
*}
	//############End Removing all cases of atherosclerosis ##############)

*}
	//########## End prepare data for analysis ###########)
	
*{	//#############  Primary Analysis(

cd `DATA_FOLDER'
use "ULSAM77 data proteomics Proseek.dta", clear 
stset allsviktdatestata, id(id) fail(allsvikt==1) enter(time beskdat70mdy) origin(time birthdate) scale(365.25)

file open table2 using "Replication_ULSAM.txt", write replace  
file write table2 "Variable" ///
		_tab "Lower 95% CI" ///
		_tab "Age Adjusted HR" ///
		_tab "upper 95% CI" ///
		_tab "P-value1"  ///
		_tab "betaP-value1" ///
		_tab "seP-value1" ///
		_tab "N_tot" ///
		_tab "N_events" ///		
	_tab "Variable" ///
		_tab "Multiple Adjusted Lower 95% CI" ///  
		_tab "Multiple Adjusted HR" ///
		_tab "Multiple Adjusted upper 95% CI" ///
		_tab "P-value2" ///
		_tab "betaP-value2" ///
		_tab "seP-value2" ///
		_tab "Multiple Adjusted N_tot" ///
		_tab "Multiple Adjusted N_events" _n

	// foreach var of varlist prot*  { // use when making supplementary table stbl01 
	foreach var of varlist prot_am prot_ccl3 prot_cd40 prot_chi3l1 prot_csf_1 prot_ctsl1 prot_cx3cl1 prot_ecp prot_fabp4 prot_fgf_23 prot_fs prot_gdf_15 prot_hgf prot_hk11 prot_il_6 prot_mcp_1 prot_mmp_10 prot_mmp_12 prot_opg prot_par_1 prot_plgf prot_scf prot_spon1 prot_st2 prot_t_pa prot_tim prot_tnf_r1 prot_tnf_r2 prot_trail_r2 prot_u_par  {

	quietly:stcox `var' 
	file write table2 "`var'" ///
		_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
		_tab (round (exp(_b[`var']),0.01)) ///
		_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
		_tab (round (2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N_sub)) ///
		_tab (e(N_fail)) ///			
		_tab
	quietly:stcox `var' exaktlder70 manuelltsbp kn medicinmothgtbt hdl ldl lipmed smoker diabetesmellitus validmi70 bmi SV2RV534 gfr afib
	file write table2 "`var'" ///
		_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
		_tab (round (exp(_b[`var']),0.01)) ///
		_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
		_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N_sub)) ///
		_tab (e(N_fail)) _n
		
}

file close table2 
insheet using "Replication_ULSAM.txt", clear
// sort multipleadjustedhr
multproc, pvalue(pvalue1) method(simes) nhcred(signpvalue1)
multproc, pvalue(pvalue2) method(simes) nhcred(signpvalue2)
g cohort="USAM"
sort variable 
// save "Stbl01 HF_ULSAM.dta", replace // This is commented out when doing primary analysis. Change "foreach var of varlist prot*" above When making supplementary table together with this line. 

keep if pvalue1 <=0.05
save "Regression table_proteomics vs HF_ULSAM.dta", replace 

*}
	//###########End Primary Analysis)

*}
	//###########End ULSAM ##############)

*{	//############# basic characteristics table and result table(
/*
Used for making result table 1 (cases of heart failure, MI, A-fib)
Table listing Basic characteristics and major cardiovascular risk factors
*/

*{	//#############  Combining bcharcs information in 4p 4u ##############(
//Table listing Basic characteristics and major heart failure risk factors
//Output file:
//"PmcsBCharcsPivusUlsam.dta"

//PIVUS 
local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
// use "PIVUS70 prevalent svikt droppade Before Stsplit.dta", clear 
use "PIVUS70 prevalent svikt droppade.dta", clear 

stset allsviktdatestata, id(id) fail(allsvikt==1) enter(time beskdat70mdy) origin(time birthdate) scale(365.25)

sort id 
g cohort="PIVUS"
g temp = id 
tostring temp, force replace
g idnew = "PIVUS"+temp
drop id 
rename idnew id 

// rename rkarenu smoker
rename medicinmothgtbt bpmed
rename validmi70 prevMI
rename SV2RV534 LVH
rename allsviktdatestata allhfdates
label variable allhfdates "Heart failure incident in Stata format"
rename allsvikt allhf
label variable allhf "Heart failure incident"
rename beskdat70mdy startdate
// rename exaktlder70 age
rename allmidatestata amidate

keep  *prot* *obnp* id cohort age bmi tg manuelltsbp manuelltdbp ldl hdl fastingglucose lipmed ntprobnp smoker bpmed LVH diabetesmellitus allhfdates allhf birthdate startdate oralaantidiabetika insulintreatment HaveAfib amidate prevMI prevAF allmi HaveMI insulintreatment oralaantidiabetika gfr male kn waist ivrt correctedef
save "ULSAM77 PIVUS70 all heart failure events.dta", replace  


//ULSAM 
use "ULSAM77 data proteomics Proseek.dta", clear 
tab insulin insulintreatment
g male = kn
label var male "Male=1"

rename medicinmothgtbt bpmed
rename validmi70 prevMI
rename SV2RV534 LVH
rename allsviktdatestata allhfdates
label variable allhfdates "Heart failure incident in Stata format"
rename allsvikt allhf
label variable allhf "Heart failure incident"
rename beskdat70mdy startdate
rename exaktlder70 age

g cohort="ULSAM"
g temp = id 
tostring temp, force replace
g idnew = "ULSAM"+temp
drop id 
rename idnew id 

sort cohort id 
// browse if HaveMI ==. // all with _t !=0 are contributing with time updated covariates 
bysort HaveMI: tab prevMI allmi // 49 individuals had MI at start of study and 94 (109 including those having prevMI after HF) at end of study, looks like more here due to stsplit

keep *prot* *obnp* id cohort age bmi tg manuelltsbp manuelltdbp ldl hdl fastingglucose lipmed ntprobnp smoker bpmed LVH diabetesmellitus allhfdates allhf birthdate startdate oralaantidiabetika insulintreatment prevAF HaveAfib amidate prevMI allmi HaveMI insulintreatment oralaantidiabetika gfr male kn waist
append using "ULSAM77 PIVUS70 all heart failure events.dta"

drop if allhfdates <= startdate

foreach var in lipmed bpmed oralaantidiabetika insulintreatment insulin diabetesmellitus { // not dropping those that have been verifiedby Johan Ärnlöv to be miscoded 
	replace `var' = 0 if `var' ==. 
	}

foreach var in kn tg hdl ldl smoker prevMI bmi gfr manuelltsbp manuelltdbp LVH lipmed fastingglucose insulintreatment oralaantidiabetika bpmed ntprobnp {
	drop if `var' ==. 
	}

g study=0 
recode study 0=1 if cohort =="ULSAM"

save "ULSAM77 PIVUS70 all heart failure events.dta", replace 
	// 2016-06-27 14.36 this is used in all secondary analysis including risk prediction model  

keep id cohort age bmi tg manuelltsbp manuelltdbp ldl hdl fastingglucose lipmed ntprobnp smoker bpmed LVH diabetesmellitus allhfdates allhf birthdate startdate oralaantidiabetika insulintreatment prevAF HaveAfib amidate prevMI allmi HaveMI insulintreatment oralaantidiabetika gfr male kn waist ivrt correctedef

save "PmcsBCharcsPivusUlsam.dta", replace 

*}
	//###########End Combining bcharcs information in 4p 4u ##############)

*{	//##############  Table 1 ##############(
 // Basic characteristics table 

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
use "PmcsBCharcsPivusUlsam.dta", clear

stset allhfdates, id(id) fail(allhf==1) enter(startdate) origin(birthdate) scale(365.25) 
duplicates drop id cohort, force // Intentionally remove time updated 
stset allhfdates, id(id) fail(allhf==1) enter(startdate) origin(birthdate) scale(365.25) 
g correctedefpercent = correctedef*100
drop id 
order age manuelltsbp manuelltdbp bmi ldl hdl tg fastingglucose ntprobnp LVH smoker lipmed bpmed male

file open table2 using "BsicTablePmcsBCharcsPivusUlsam.txt", replace write
file write table2 "Cohort (Number of individuals)" _tab
	foreach cohort in PIVUS ULSAM {
		file write table2 "`cohort'"
		sum age if cohort=="`cohort'"
		file write table2 " (" (r(N)) ")" _tab 
		}		
	file write table2 _n 
	foreach var of varlist age {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS ULSAM {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),0.1)) " (" (round(r(p25),0.1)) "-" (round(r(p75),0.1)) ")"  _tab
			}
		file write table2 _n
		}

	foreach var of varlist bmi {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS ULSAM {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),0.1)) " (" (round(r(p25),0.1)) "-" (round(r(p75),0.1)) ")"  _tab
			} 
		file write table2 _n
		}

	foreach var of varlist ldl hdl tg fastingglucose {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS ULSAM {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),0.1)) " (" (round(r(p25),0.1)) "-" (round(r(p75),0.1)) ")"  _tab
			}
		file write table2 _n
		}

	foreach var of varlist manuelltsbp manuelltdbp gfr {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS ULSAM {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),1)) " (" (round(r(p25),1)) "-" (round(r(p75),1)) ")"  _tab
			}
		file write table2 _n
		}

	foreach var of varlist ntprobnp ivrt {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),1)) " (" (round(r(p25),1)) "-" (round(r(p75),1)) ")"  _tab
			}
		foreach cohort in ULSAM {
			file write table2 "NA" 
			}
		file write table2 _n
		}

	foreach var of varlist correctedefpercent {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS {
			sum `var' if cohort=="`cohort'", d
			file write table2 (round(r(p50),1)) " (" (round(r(p25),1)) "-" (round(r(p75),1)) ")"  _tab
			}
		foreach cohort in ULSAM {
			file write table2 "NA" 
			}
		file write table2 _n
		}

	foreach var of varlist male smoker LVH lipmed bpmed prevMI prevAF diabetesmellitus insulintreatment oralaantidiabetika {
		file write table2 "`var'" _tab
		foreach cohort in PIVUS ULSAM {
			sum `var' if cohort=="`cohort'" 
				local x=r(N) 
			sum `var' if cohort=="`cohort'" & `var' ==1
			/* only needed when you have >2 categories 
						// local y=r(N) 
					// sum `var' if cohort=="`cohort'" & `var' ==2
						// local z=r(N) 
			*/
			file write table2 (r(N)) " (" (round(r(N)/`x'*100),1) "%)" _tab
			}
		file write table2 _n
		}

file close table2  
insheet using "BsicTablePmcsBCharcsPivusUlsam.txt", clear

foreach var in v2 {
	tostring `var', force replace format("%04.2f")
		// to get it in e.g. the following format "0.12"
	}

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("T1")

*}
	//############End Table 1 ############# )
	
*}
	//##########End basic characteristics table and result tables)

*{	//#############  Forest plot ##############(
// 2016-02-09 14.22 
	// tried making code below more efficient, but the current method appears to be esiest to expand on later so keeping it as it is
// 2016-01-27 09.40 
	// updating code below for new unadjusted analysis 
	// structured the Cox-regression so that it is easier to switch between the two analysis methods 
// 2015-10-21 18.01 
	// this is best to use for 3ua01
	// Latest general version is in template Stata Do-file r726-800 2015-10-21 10.24 

*{	//##############  fplot significant proteins ##############(
local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'

//ULSAM 
import excel using "c:/Dropbox/Work/Uppsala/Various/OlinkInfo/ConvertionTable/primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 

merge 1:1 variable using "Regression table_proteomics vs HF_ULSAM.dta"
keep if _merge==3
drop _merge
sort variable 

gen cirange=string(lower95ci)+", "+string(upper95ci)
keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 cohort cirange
order name_and_abbreviation variable UniProtID name abbreviation lower95ci ageadjustedhr upper95ci pvalue1 cohort cirange
save "Regression table_proteomics vs svikt_PIVUS_ULSAM.dta", replace
gsort -ageadjustedhr variable

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("f14u")

order abbreviation lower95ci ageadjustedhr upper95ci pvalue1
keep abbreviation lower95ci ageadjustedhr upper95ci pvalue1

export excel using "c:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/PublicationResults/STblFplotTop18ULSAM/STblFplotTop18ULSAM.xlsx", sheetreplace firstrow(var) sheet("STblFplotTop18ULSAM") 


//PIVUS 
cd `DATA_FOLDER'
import excel using "c:/Dropbox/Work/Uppsala/Various/OlinkInfo/ConvertionTable/primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
sort variable 
merge variable using "Regression table_proteomics vs svikt_PIVUS_ULSAM.dta" // NB! - merging ULSAM and PIVUS to get same sort 
keep if _merge==3
drop _merge
sort variable 
g ULSAMHR = ageadjustedhr
	// in order to get identical sort in both cohorts 
keep name_and_abbreviation abbreviation variable UniProtID name cohort ULSAMHR
merge variable using "Regression table_proteomics vs svikt_PIVUS_HR.dta"
keep if _merge==3
drop _merge
gsort -ULSAMHR variable 
keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 
g cohort="PIVUS"
order name_and_abbreviation variable UniProtID name abbreviation lower95ci ageadjustedhr upper95ci pvalue1 cohort 
gen cirange=string(lower95ci)+", "+string(upper95ci)

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("f14p") 

order abbreviation lower95ci ageadjustedhr upper95ci pvalue1
keep abbreviation lower95ci ageadjustedhr upper95ci pvalue1

export excel using "c:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/PublicationResults/STblFplotTop18PIVUS/STblFplotTop18PIVUS.xlsx", sheetreplace firstrow(var) sheet("STblFplotTop18PIVUS")

*}
	//############End fplot significant proteins ##############)

*{	//#############  forest plot on all proteins ##############(
//ULSAM 
cd `DATA_FOLDER'
import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 

merge 1:1 variable using "Stbl01 HF_ULSAM.dta"
keep if _merge==3
drop _merge
gen cirange=string(lower95ci)+", "+string(upper95ci)
keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 betapvalue1 sepvalue1 multipleadjustedlower95ci multipleadjustedhr multipleadjustedupper95ci cirange cohort
order name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 betapvalue1 sepvalue1 multipleadjustedlower95ci multipleadjustedhr multipleadjustedupper95ci cirange cohort
sort variable 
save "Stbl01 HF_PIVUS_ULSAM.dta", replace
gsort -ageadjustedhr variable
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("4uAllPrts")

order abbreviation lower95ci ageadjustedhr upper95ci pvalue1
keep abbreviation lower95ci ageadjustedhr upper95ci pvalue1
export excel using "c:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/PublicationResults/STblFplotAllULSAM/STblFplotAllULSAM.xlsx", sheetreplace firstrow(var) sheet("STblFplotAllULSAM") 

//PIVUS 
cd `DATA_FOLDER'
import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
sort variable 
merge variable using "Stbl01 HF_PIVUS_ULSAM.dta" // Merging ULSAM and PIVUS to get same sort 
keep if _merge==3
drop _merge
sort variable 

g ULSAMHR = ageadjustedhr // in order to get identical sort in both cohorts 
keep name_and_abbreviation abbreviation variable UniProtID name cohort ULSAMHR
merge variable using "Regression table_proteomics vs svikt_PIVUS_HR.dta"
keep if _merge==3
drop _merge
gsort -ULSAMHR variable 

gen cirange=string(lower95ci)+", "+string(upper95ci)
keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 betapvalue1 sepvalue1 multipleadjustedlower95ci multipleadjustedhr multipleadjustedupper95ci cirange 
order name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedhr upper95ci pvalue1 betapvalue1 sepvalue1 multipleadjustedlower95ci multipleadjustedhr multipleadjustedupper95ci cirange 
g cohort="PIVUS"

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("4pAllPrts") 

order abbreviation lower95ci ageadjustedhr upper95ci pvalue1
keep abbreviation lower95ci ageadjustedhr upper95ci pvalue1
export excel using "c:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/PublicationResults/STblFplotAllPIVUS/STblFplotAllPIVUS.xlsx", sheetreplace firstrow(var) sheet("STblFplotAllPIVUS") 

*}
	//###########End forest plot on all proteins ##############)

*}
	//###########End Forest plot ##############)

*{ 	//############# meta-analysis ##############(
/*
2016-06-08 12.05 
	removed NT-proBNP (alone and with fully adjusted)
2016-02-16 17.02 
	has time dependent covariates "HaveAfib" and "HaveMI" 
*/ 

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
use "ULSAM77 PIVUS70 all heart failure events.dta", clear 

foreach var of varlist male bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp bpmed smoker LVH gfr HaveMI HaveAfib ntprobnp {
	drop if `var' ==. 
	}

stset allhfdates, id(id) fail(allhf==1) enter(time startdate) origin(time birthdate) scale(365.25) 

file open table2 using "meta-analysis.txt", write replace  
file write table2 "Variable" ///
		_tab "Lower 95% CI" ///
		_tab "Age Adjusted HR" ///
		_tab "upper 95% CI" ///
		_tab "P-value1"  ///
		_tab "betaP-value1" ///
		_tab "seP-value1" ///
		_tab "N_tot" ///
		_tab "N_events" ///		
		_tab "Protein Counter" ///		
		_tab "Model Counter" ///		
	_tab "Model" _n ///

local j=0 

	foreach var of varlist prot_gdf_15 prot_am prot_fgf_23 prot_tim prot_trail_r2 prot_spon1 prot_mmp_12 prot_fs prot_u_par prot_plgf prot_fabp4 prot_csf_1 prot_tnf_r1 prot_opg prot_st2 prot_par_1 prot_chi3l1 prot_tnf_r2 {

local i=1 
local j=`j'+1 
		quietly:stcox `var' male study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///			
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Baseline model" /// 
			_n

				local i=`i'+1 
		quietly:stcox `var' male bmi study
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "BMI" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male ldl hdl tg lipmed study
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Lipids" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male fastingglucose insulintreatment oralaantidiabetika study
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Glycemic status" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male manuelltsbp manuelltdbp bpmed study
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Blood pressure" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male smoker study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Smoking" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male LVH study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Left ventricular hypertrophy" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male HaveMI study  // using time updated, not prevMI prevAF
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Myocardial infarction prior/during study" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male HaveAfib study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Atrial fibrillation prior/during study" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male gfr study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Kidney function" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male stdntprobnp study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "NT-proBNP (PIVUS ELISA ULSAM Olink)" ///
			_n

				local i=`i'+1 
		quietly:stcox `var' male bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp bpmed smoker LVH gfr HaveAfib HaveMI study
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Fully adjusted" ///
			_n
 
				local i=`i'+1 
		quietly:stcox `var' male bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp bpmed smoker LVH gfr HaveMI HaveAfib stdntprobnp study 
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "`i'" /// 
			_tab "Fully adjusted + NT-proBNP (PIVUS ELISA ULSAM Olink)" ///
			_n
		}

file close table2
insheet using "meta-analysis.txt", clear
multproc, pvalue(pvalue1) method(simes) nhcred(signpvalue1)
save "Results.dta", replace 

import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
merge m:m variable using "Results.dta"

keep if _merge==3
drop _merge

sort proteincounter modelcounter 
order name_and_abbreviation variable UniProtID name abbreviation lower95ci ageadjustedhr upper95ci pvalue1
gen cirange=string(lower95ci)+", "+string(upper95ci)

save "Results.dta", replace 

drop if model =="Fully adjusted + NT-proBNP (PIVUS ELISA ULSAM Olink)" 
drop if model =="NT-proBNP (PIVUS ELISA ULSAM Olink)" 

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("heatmap") 

use "Results.dta", clear
sort proteincounter modelcounter 
order abbreviation model lower95ci ageadjustedhr upper95ci pvalue1
keep abbreviation model lower95ci ageadjustedhr upper95ci pvalue1

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/PublicationResults/STblHeatmap/Supplemental Table heatmap.xlsx", sheetreplace firstrow(var) sheet("STblHeatmap") 

use "Results.dta", clear
keep if model =="Baseline model" 
rename ageadjustedhr BaselineModel 
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("BaselineModel") 

use "Results.dta", clear
keep if model =="Fully adjusted" 
rename ageadjustedhr FullyAdjusted
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("FullyAdjusted") 

use "Results.dta", clear
keep if model =="Fully adjusted + NT-proBNP (PIVUS ELISA ULSAM Olink)" 
keep if inlist(variable, "prot_fs", "prot_gdf_15", "prot_mmp_12", "prot_opg", "prot_par_1", "prot_st2", "prot_tim", "prot_trail_r2", "prot_u_par")
rename ageadjustedhr FullyAdjustedNTproBNP
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("FullyAdjustedNTproBNP") 

*}
	//######### END meta-analysis ##############)

*{	//#############  LVFunction ##############(

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
// use "PIVUS70 prevalent svikt droppade.dta", clear // has time updated covariates 
use "PIVUS70 prevalent svikt droppade Before Stsplit.dta", clear

stset allsviktdatestata, fail(allsvikt) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)

tabstat ivrt, stat (mean sd p50) // p50 is median 	
tabstat correctedef, stat (mean sd p50) 
sort correctedef

file open table2 using "Table 2 LV function.txt", write replace  
file write table2 "Protein name" ///
		_tab "CI" /// 
		_tab "Lower 95% CI" ///
		_tab "Age Adjusted beta Systolic" ///
		_tab "Upper 95% CI" ///
		_tab "P-value1 Systolic" ///
		_tab "-value1" ///
		_tab "seP-value1" ///
		_tab "N_tot" ///
		_tab "ProteinOrder" ///
		_tab "NamesAgain" ///		
	_tab "CIDiastolic" ///
		_tab "Lower 95% CI Diastolic" ///
		_tab "Age Adjusted beta Diastolic" ///
		_tab "Upper 95% CI Diastolic" ///
		_tab "P-value1 Diastolic" ///
		_tab "betaP-value1 Diastolic" ///
		_tab "seP-value1 Diastolic" ///
		_tab "N_tot Diastolic" _n

local j=0	

foreach var of varlist prot_gdf_15 prot_tim prot_trail_r2 prot_spon1 prot_mmp_12 prot_fs prot_u_par prot_opg prot_st2 {

local j=`j'+1 
	
		// _tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		// _tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
	regress correctedef `var' exaktlder70 kn
	file write table2 "`var'" ///
		_tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		_tab (round (100*_b[`var']-100*1.96*_se[`var']),0.001) ///
		_tab (round (100*_b[`var'],0.001)) ///
		_tab (round (100*_b[`var']+100*1.96*_se[`var']),0.001) ///
		_tab (round (2*(ttail(e(df_r),abs(_b[`var']/_se[`var']))),0.0000000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N)) ///
		_tab "`j'" ///
		_tab 

	regress ivrt `var' exaktlder70 kn
	file write table2 "`var'" ///
		_tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		_tab (round (100*_b[`var']-100*1.96*_se[`var']),0.001) ///
		_tab (round (100*_b[`var'],0.001)) ///
		_tab (round (100*_b[`var']+100*1.96*_se[`var']),0.001) ///
		_tab (round (2*(ttail(e(df_r),abs(_b[`var']/_se[`var']))),0.0000000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N)) ///
		_n
	
	}

file close table2 
insheet using "Table 2 LV function.txt", clear 
rename proteinname variable 

save "Results.dta", replace 
import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
merge m:m variable using "Results.dta"

keep if _merge==3
drop _merge

gen cirange=string(lower95ci)+", "+string(upper95ci)
// keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedbetasystolic upper95ci pvalue1systolic cirange
order name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedbetasystolic upper95ci pvalue1systolic cirange

sort proteinorder
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx",  sheetreplace firstrow(var) sheet("systolicDiastolic") 

*}
	//###########End LVFunction ##############)

*{ 	//############# Sensitivity analysis of IVRT ##############(
/*
On MMP-12 and TIM-1: Drop all EF <0.5 (to make sure that we don't pseudo normalize) and perform Cox-regression on this ER quote 
2016-09-22 09.45	updated analysis 
*/ 

local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
// use "PIVUS70 prevalent svikt droppade.dta", clear // has time updated covariates 
use "PIVUS70 prevalent svikt droppade Before Stsplit.dta", clear

stset allsviktdatestata, fail(allsvikt) id(id) origin(birthdate) enter(beskdat70mdy) scale(365.25)
// duplicates drop id cohort, force 

foreach var of varlist male bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp smoker gfr ntprobnp {
	drop if `var' ==. 
	}

drop if correctedef <.5 

file open table2 using "Table 2 LV function.txt", write replace  
file write table2 "Protein name" ///
		_tab "CI" /// 
		_tab "Lower 95% CI" ///
		_tab "Age Adjusted beta Systolic" ///
		_tab "Upper 95% CI" ///
		_tab "P-value1 Systolic" ///
		_tab "-value1" ///
		_tab "seP-value1" ///
		_tab "N_tot" ///
		_tab "ProteinOrder" ///
		_tab "NamesAgain" ///		
	_tab "CIDiastolic" ///
		_tab "Lower 95% CI Diastolic" ///
		_tab "Age Adjusted beta Diastolic" ///
		_tab "Upper 95% CI Diastolic" ///
		_tab "P-value1 Diastolic" ///
		_tab "betaP-value1 Diastolic" ///
		_tab "seP-value1 Diastolic" ///
		_tab "N_tot Diastolic" _n

local j=0	
foreach var of varlist prot_tim prot_mmp_12 {

local j=`j'+1 
	
		// _tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		// _tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
	regress correctedef `var' exaktlder70 kn
	file write table2 "`var'" ///
		_tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		_tab (round (100*_b[`var']-100*1.96*_se[`var']),0.001) ///
		_tab (round (100*_b[`var'],0.001)) ///
		_tab (round (100*_b[`var']+100*1.96*_se[`var']),0.001) ///
		_tab (round (2*(ttail(e(df_r),abs(_b[`var']/_se[`var']))),0.0000000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N)) ///
		_tab "`j'" ///
		_tab 

	regress ivrt `var' exaktlder70 kn
	file write table2 "`var'" ///
		_tab " (" (round (_b[`var']-1.96*_se[`var']),0.001) ", "  (round (_b[`var']+1.96*_se[`var']),0.001) ")"  ///
		_tab (round (100*_b[`var']-100*1.96*_se[`var']),0.001) ///
		_tab (round (100*_b[`var'],0.001)) ///
		_tab (round (100*_b[`var']+100*1.96*_se[`var']),0.001) ///
		_tab (round (2*(ttail(e(df_r),abs(_b[`var']/_se[`var']))),0.0000000001)) ///
		_tab (_b[`var']) ///
		_tab (_se[`var']) ///
		_tab (e(N)) ///
		_n
	
	}

file close table2 
insheet using "Table 2 LV function.txt", clear 
rename proteinname variable 

save "Results.dta", replace 
import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
merge m:m variable using "Results.dta"

keep if _merge==3
drop _merge

gen cirange=string(lower95ci)+", "+string(upper95ci)
// keep name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedbetasystolic upper95ci pvalue1systolic cirange
order name_and_abbreviation abbreviation variable UniProtID name lower95ci ageadjustedbetasystolic upper95ci pvalue1systolic cirange

sort proteinorder
export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx",  sheetreplace firstrow(var) sheet("SysDiasEFB0.5") 

*}
	//######### END Sensitivity analysis of IVRT ##############)

*{ 	//############# Sensitivity analysis of NIHF ##############(
/*
nonischemic heart failure (NIHF) sensitivity analysis 
Perform stratified analysis on top 9 
*/ 

*{	//#############  merge PIVUS and ULSAM prior to stsplit ##############(
//PIVUS 
local DATA_FOLDER = "c:\Dropbox\Work\Uppsala\0ToDo\08NotSynced\Pr1\Pr1Datafiles\"
cd `DATA_FOLDER'
use "PIVUS70 prevalent svikt droppade Before Stsplit.dta", clear 

stset allsviktdatestata, id(id) fail(allsvikt==1) enter(time beskdat70mdy) origin(time birthdate) scale(365.25)

sort id 
g cohort="PIVUS"
g temp = id 
tostring temp, force replace
g idnew = "PIVUS"+temp
drop id 
rename idnew id 

// rename rkarenu smoker
rename medicinmothgtbt bpmed
rename validmi70 prevMI
rename SV2RV534 LVH
rename allsviktdatestata allhfdates
label variable allhfdates "Heart failure incident in Stata format"
rename allsvikt allhf
label variable allhf "Heart failure incident"
rename beskdat70mdy startdate
// rename exaktlder70 age
rename allmidatestata amidate

keep  *prot* *obnp* id cohort age bmi tg manuelltsbp manuelltdbp ldl hdl fastingglucose lipmed ntprobnp smoker bpmed LVH diabetesmellitus allhfdates allhf birthdate startdate oralaantidiabetika insulintreatment amidate prevMI prevAF allmi insulintreatment oralaantidiabetika gfr male kn waist ivrt correctedef
save "PIVUS70 all heart failure events Before stsplit temporary.dta", replace  

//ULSAM 
use "ULSAM77 data proteomics Proseek Before Stsplit.dta", clear 
tab insulin insulintreatment
g male = kn
label var male "Male=1"

rename medicinmothgtbt bpmed
rename validmi70 prevMI
rename SV2RV534 LVH
rename allsviktdatestata allhfdates
label variable allhfdates "Heart failure incident in Stata format"
rename allsvikt allhf
label variable allhf "Heart failure incident"
rename beskdat70mdy startdate
rename exaktlder70 age

g cohort="ULSAM"
g temp = id 
tostring temp, force replace
g idnew = "ULSAM"+temp
drop id 
rename idnew id 

sort cohort id 

keep *prot* *obnp* id cohort age bmi tg manuelltsbp manuelltdbp ldl hdl fastingglucose lipmed ntprobnp smoker bpmed LVH diabetesmellitus allhfdates allhf birthdate startdate oralaantidiabetika insulintreatment prevAF amidate prevMI allmi  insulintreatment oralaantidiabetika gfr male kn waist
append using "PIVUS70 all heart failure events Before stsplit temporary.dta"

g study=0 
recode study 0=1 if cohort =="ULSAM"

save "ULSAM77 PIVUS70 all heart failure events Before stsplit.dta", replace 
*}
	//#########  end merge PIVUS and ULSAM prior to stsplit ##############)

// use "ULSAM77 PIVUS70 all heart failure events.dta", clear 
// use "ULSAM77 PIVUS70 all heart failure events Before stsplit.dta", clear 
sort allhf amidate 

drop if prevMI == 1 // remove 121 individuals with validmi70 at or before start of study 
drop if amidate < startdate
// browse id prevMI startdate amidate allhfdates allhf birthdate if amidate < startdate
// tab allmi if amidate < startdate

tab allmi allhf 
sort cohort allhfdates
tab allmi if amidate != . 
recode allhf 1=0 if amidate <= allhfdates & allmi == 1 // 24 changes 

replace allhfdates = amidate if amidate < allhfdates & allmi == 1 // 75 changes (some individuals with non-validmi70  prior to study will have outcome before start of study)
stset allhfdates, fail(allhf) id(id) origin(birthdate) enter(startdate) scale(365.25) 

file open table2 using "Results.txt", write replace  
file write table2 "Variable" ///
		_tab "Lower 95% CI" ///
		_tab "HR" ///
		_tab "upper 95% CI" ///
		_tab "P-value1"  ///
		_tab "beta" ///
		_tab "se" ///
		_tab "N_tot" ///
		_tab "N_events" ///		
		_tab "Protein Counter" ///		
	_tab "Model" _n ///

local j=0	
	foreach var of varlist prot_fs prot_gdf_15 prot_mmp_12 prot_opg prot_par_1 prot_st2 prot_tim prot_trail_r2 prot_u_par {

local j=`j'+1 
		quietly:stcox `var' male bmi ldl hdl tg lipmed fastingglucose insulintreatment oralaantidiabetika manuelltsbp manuelltdbp bpmed smoker LVH gfr study prevAF
		file write table2 "`var'" ///
			_tab (round (exp(_b[`var']-1.96*_se[`var']),0.01)) ///
			_tab (round (exp(_b[`var']),0.01)) ///
			_tab (round (exp(_b[`var']+1.96*_se[`var']),0.01)) ///
			_tab (round(2*(1-normal(abs(_b[`var']/_se[`var']))),0.00000001)) ///
			_tab (_b[`var']) ///
			_tab (_se[`var']) ///
			_tab (e(N_sub)) ///
			_tab (e(N_fail)) ///
			_tab "`j'" /// 
			_tab "Fully adjusted" ///
			_n
		}

file close table2
insheet using "Results.txt", clear
multproc, pvalue(pvalue1) method(simes) nhcred(signpvalue1) 
save "Results.dta", replace 

import excel using "c:\Dropbox\Work\Uppsala\Various\OlinkInfo\ConvertionTable\primary_key_olink_variable_name_full_name_UniProt.xlsx", firstrow clear 
merge m:m variable using "Results.dta"

keep if _merge==3
drop _merge

sort hr
order name_and_abbreviation variable UniProtID name abbreviation lower95ci hr upper95ci pvalue1
gen cirange=string(lower95ci)+", "+string(upper95ci)

save "Results.dta", replace 

export excel using "C:/Dropbox/Work/Uppsala/0ToDo/project01/ResultsFiles/Main results/results 3up01 pek01resus 3ua01 ULSAM PIVUS cox-regression multproc main results mresus.xlsx", sheetreplace firstrow(var) sheet("NIHF") 

*}
	//######### END Sensitivity analysis of NIHF ##############)

*{	//#############  reference information ##############(

*{	//##############  study information ##############(
// proteomicsHF version 2016-12-09 
// primary analysis 
	//discovery in PIVUS 
		//Cox proportional hazards model adjusting for age and sex 
		//5 % FDR multiple testing cutoff 
	//Replicate in ULSAM 
		//nominal p-value (0.05) cutoff 

// Secondary analyses 
	// combined cohort ("Meta-analysis") and adjust for each confounder 
		// Visualization using heatmap  
	// Proteins that remain significant in fully adjustedadjusted 
		// forest plot 
	// Determine which biomarkers are associated to systolic and/or diastolic function 
	// Prediction model
		// separate R document with Atherosclerosis Risk in Communities (ARIC) score 
		// perform for 
			// all proteins 

*}
	//############End study information ##############)

*}
	//########## End reference information ###########)