##########################################################################
#################### Step 1: data import and cleanup  ####################
####################        Madeline and Carina      #####################
##########################################################################


####LIBRARIES AND FILE PATHS####
source('T:\\Analysis\\Madeline\\MS_Precip_Libraries.R')

#Load if existing
load(paste0(f1,'GI_Resp_Anal.RData'))

####IMPORT ANCILLARY FILES####
#ZIP to ZCTA crosswalk 
zip_to_zcta_mi_oh_pa <- read_excel(paste0(f1,"zip_to_zcta_mi_oh_pa.xlsx"))
zip_to_zcta_mi_oh_pa$ZIP_CD<-as.integer(zip_to_zcta_mi_oh_pa$ZIP_CD)
zip_to_zcta_mi_oh_pa$ZCTA<-as.integer(zip_to_zcta_mi_oh_pa$ZCTA)

####IMPORTING CLAIMS DATA AND SETTING TO DT####
#Note: Some lines are commented out because we need to protect the files from being overwritten or because we don't use those lines on every year of data.
#The following criteria were used to subset the original data and were also applied to the denominator.
#AGE_AT_END_REF_YR = 65 or older
#STATE_CNTY_FIPS_CD_{x} = a code in the researcher submitted file (must meet this criteria for at least one month, does not need to meet in all months of the year)
#HMO_IND_{x} = ("0","1","4","A") for all months of the year OR until death OR for all months from first enrollment month of the year through the end of the year or until death (for example, if a bene enrolls in September 2016, must meet this criteria for September-December 2016 or until death to qualify)
#MDCR_ENTLMT_BUYIN_IND_{x} = 3 or C for all months of the year OR until death OR for all months from first enrollment month of the year through the end of the year or until death (for example, if a bene enrolls in September 2016, must meet this criteria for September-December 2016 or until death to qualify)

#claims_dim<-data.table(year=2006:2013, gi.resp.subset.import=as.numeric(NA), mx.race.cd.per.bene_id=as.numeric(NA), mx.sex.cd.per.bene_id=as.numeric(NA), mx.gndr.cd.per.bene_id=as.numeric(NA), nodup.present=as.numeric(NA), age65up=as.numeric(NA), zip.to.zcta.present=as.numeric(NA),prism.present=as.numeric(NA))
#claims_master<-data.table()


#vectors used in the subsequent loop
icdCols<-c("PRNCPAL_DGNS_CD",paste0('ICD_DGNS_CD',1:11),'ICD_DGNS_E_CD12',paste0('DGNS_',1:10,'_CD'))
claimsCols<-c("BENE_ID","ZIP_CD", "BENE_BIRTH_DT","DOB_DT",'COVSTART',"CLM_FROM_DT",'CLM_THRU_DT',"AGE_AT_END_REF_YR", "GNDR_CD", "SEX_IDENT_CD", "BENE_RACE_CD", "RTI_RACE_CD",'SOURCE', "CLM_LINE_NUM", "BENE_SEX_CD", icdCols, "gi_disease","resp_disease")
viralresp_cd<-c("460", "461", "462", "463", "464", "465", "466", "480", "487", "488")
viralgi_cd<-c("008", "009")

#identify the column names in each year and see if any are missing
#colNamesSAS<-list()
for (yy in 2006:2013) {
  cat(yy)
  v1<-paste0('claims_gi_resp_sub_insample_',yy,'.sas7bdat')
  claims1<-read_sas(paste0(f2,'\\',v1),n=1)
  colNamesSAS[[as.character(yy)]]<-names(claims1)
  cat(' missing columns:',claimsCols[!(claimsCols %in% names(claims1))],'\n')
}

for (yy in 2006:2013) {print(claimsCols[!(claimsCols %in% colNamesSAS[[as.character(yy)]])])}
#none missing

for (yy in 2006:2013) {
cat('\n--',yy,'--\n')
v1<-paste0('claims_gi_resp_sub_insample_',yy,'.sas7bdat')
v3<-claimsCols[claimsCols %in% colNamesSAS[[as.character(yy)]]]
claims1<-read_sas(paste0(f2,'\\',v1),col_select=all_of(v3))
rm(v1)

#Converting to data table
claims1<-setDT(claims1)

#Record initial number of rows
claims_dim[year==yy,gi.resp.subset.import:=nrow(claims1)]

####Does BENE_ID ever have multiple race, SEX, or BIRTH_DTs?####
v1<-paste0(claims1$BENE_ID,claims1$BENE_RACE_CD)
claims_dim[year==yy,mx.race.cd.per.bene_id:=length(unique(v1))-length(unique(claims1$BENE_ID))]
#0 = Unknown, 1 = White, 2 = Black, 3 = Other, 4 = Asian, 5 = Hispanic, 6 = North American Native
#All missing in 2006 data, so no race analyses possible with the 2006 data unless the individual appears in subsequent years.
#A lot of non-matching race codes in the 2007 2013 data. Look at 2013:
#v5<-claims1$BENE_ID[duplicated(claims1$BENE_ID)]
#v3<-v1[!duplicated(v1)]
#v4<-claims1$BENE_ID[claims1$BENE_ID %in% v5 & v1 %in% v3]
#claims1[claims1$BENE_ID %in% v4[1:10],] #all instances where race code was missing for one of the entries
#length(unique(v1[claims1$BENE_RACE_CD!='']))-length(unique(claims1$BENE_ID[claims1$BENE_RACE_CD!='']))
#Only 126 instances and 8 instances of non-matching race codes where one wasn't simply missing in the 2007 and 2013 data, respectively.
#If using race in the analyses, make sure to first aggregate by BENE_ID and race, dropping missings, and then re-assign race through a merge with the aggregated data set so as to ensure that race information is retained from rows that might otherwise have been dropped due to duplicate disease types on the same day. This might also pick up race information to apply to the 2006 data.
#rm(v1,v5,v3,v4)
rm(v3)

v1<-paste0(claims1$BENE_ID,claims1$BENE_SEX_CD)
claims_dim[year==yy,mx.sex.cd.per.bene_id:=length(unique(v1))-length(unique(claims1$BENE_ID))]

v1<-paste0(claims1$BENE_ID,claims1$GNDR_CD)
claims_dim[year==yy,mx.gndr.cd.per.bene_id:=length(unique(v1))-length(unique(claims1$BENE_ID))]

v1<-paste0(claims1$BENE_ID,claims1$SEX_IDENT_CD)
claims_dim[year==yy,mx.sexident.cd.per.bene_id:=length(unique(v1))-length(unique(claims1$BENE_ID))]
claims1[,table(SEX_IDENT_CD,useNA='a')]
#Use GNDR instead of SEX for all but 2007.
#Many GNDR missings in 2006-2008 data, too. Use SEX_IDENT_CD instead in these years.

#In 2006, 6% of the SEX/GENDER identities are inconsistent within BENE_ID. Considering that Robert had to crosswalk BENE_IDs from the MedPAR records to the outpatient records, are these inconsistencies primarily in the MedPAR records?
if (yy==2006) {
  c1<-claims1[,SOURCE=='medpar2006_rev']
  v1<-paste0(claims1$BENE_ID[c1],claims1$GNDR_CD[c1])
  length(unique(v1))-length(unique(claims1$BENE_ID[c1]))
  #Only 2 of the records have GNDR_CD mismatches in the MedPAR records.

  #Is the explanation for the 21.6% GNDR_CD mismatch by BENE_ID due to a different BENE_ID system being used for the two outpatient visit requests (DUA 51193 and DUA 52573)?
  c1<-claims1[,SOURCE=='outpatient_revenue_cente']
  v1<-paste0(claims1$BENE_ID[c1],claims1$CLM_LINE_NUM[c1])
  length(unique(v1))-length(unique(claims1$BENE_ID[c1]))
  #No.

  #Is BENE_SEX_CD (from MBSF) consistent within BENE_IDs?
  c1<-claims1[,SOURCE=='outpatient_revenue_cente']
  v1<-paste0(claims1$BENE_ID[c1],claims1$BENE_SEX_CD[c1])
  length(unique(v1))-length(unique(claims1$BENE_ID[c1]))
  #Yes. There are 0 mismatches within BENE_IDs for the BENE_SEX_CD, which comes from MBSF.

  #Is DOB consistent between claims and MBSF?
  c1<-claims1[,SOURCE=='outpatient_revenue_cente']
  v1<-paste0(claims1$BENE_ID[c1],claims1$BENE_BIRTH_DT[c1])
  length(unique(v1))-length(unique(claims1$BENE_ID[c1]))
  #0
  v1<-paste0(claims1$BENE_ID[!c1],claims1$DOB_DT[!c1])
  length(unique(v1))-length(unique(claims1$BENE_ID[!c1]))                          
  #0, So DOB is consisten between claims and MBSF.
  rm(c1,v1)
}

####DROPPING DUPLICATES####
#It's not uncommon for there to be more than one claim row for a single visit date. This code drops instances where the disease type listed in the claim rows is the same.
claims1<-setorder(claims1,BENE_ID,CLM_FROM_DT,CLM_THRU_DT,CLM_LINE_NUM)
nodup_claims1 <-claims1[!duplicated(claims1[, c('BENE_ID','CLM_FROM_DT','CLM_THRU_DT','resp_disease','gi_disease')]), ]
claims_dim[year==yy,nodup.present:=nrow(nodup_claims1)]
 
#2009: 88.25% of the claims remain.
rm(claims1)

#does this make sense in 2009?
#table(claims1$CLM_LINE_NUM,useNA='a')
#sum(is.na(claims1$CLM_LINE_NUM) | claims1$CLM_LINE_NUM==1)/nrow(claims1)*100
#71% of the records do not have multiple claims for a single visit, so since this number is smaller than the 88.25% above, we can have some confidence that we are not losing too many claims.

#CHECKING FOR NAS
#cat('instances where both gi_disease and resp_disease are NA:',nrow(nodup_claims1[which(is.na(nodup_claims1$gi_disease)& (is.na(nodup_claims1$resp_disease))),]),'\n')

####CHECKING IF CLAIMS HAVE ANYONE W/ AGE<65 ####
#Three steps: 
#Step 1 checking dates are dates and converting if not 
#Step 2 creating new age variable called NEW_AGE 
#Step 3 dropping people where NEW_AGE<65

#BENE_AGE_CNT is the beneficiary's age as of date of admission
#Creating new age variable 
nodup_claims1$NEW_AGE<-as.integer(round(((as.Date(nodup_claims1$CLM_FROM_DT)-as.Date(nodup_claims1$BENE_BIRTH_DT))/365),0))
#summary(as.numeric(nodup_claims1$NEW_AGE))

#drop the individuals who weren't 65 by the time of their date of admission
nodup_claims1<-nodup_claims1[NEW_AGE>=65,]
claims_dim[year==yy,age65up:=nrow(nodup_claims1)]

#####IDENTIFYING VIRAL AND RESP CAUSES#####
v1<-apply(nodup_claims1[,..icdCols],1,function(x) any( substr(x,1,3) %in% viralresp_cd ))
nodup_claims1[,viralresp:=v1]

v3<-apply(nodup_claims1[,..icdCols],1,function(x) any( substr(x,1,3) %in% viralgi_cd ))
nodup_claims1[,viralgi:=v3]

rm(v1,v3)

####IDENTIFY GI PRINCIPAL AND RESPIRATORY PRINCIPAL####


####DROPPING COLUMNS FROM CLAIMS DATA####
#Using the CLM_FROM_DT variable for date not ADMSN_DT--ADMSN_DT has NAs
claims_red<-nodup_claims1[,c("BENE_ID", 'ZIP_CD', "CLM_FROM_DT", "CLM_THRU_DT","SEX_IDENT_CD","NEW_AGE","BENE_RACE_CD","gi_disease", "resp_disease", "viralresp", "viralgi")]

#convert numeric to integer values to save space
numCols<-c('SEX_IDENT_CD','NEW_AGE','BENE_RACE_CD','gi_disease','resp_disease','viralresp','viralgi')
claims_red<-claims_red[,(numCols):= lapply(.SD,as.integer), .SDcols=numCols]
rm(numCols)

####AGGREGATE DATA AND DROP PRE-2006 DATA####
#The annual files are based on CLM_TO_DATE but we aggregated by CLM_FROM_DATE, so there might be duplicates in the stacked master file. We also drop dates occurring before 2006-01-01.
#If we want to look at effect modification by age, race or sex, these groupings will also need to be created and added here.
claims_agg<-claims_red[CLM_FROM_DT>=as.Date('2006-01-01'), lapply(.SD, sum, na.rm = TRUE), by = c("ZIP_CD", "CLM_FROM_DT"), .SDcols =c("gi_disease", "resp_disease", "viralresp", "viralgi")]
numCols<-c('ZIP_CD','gi_disease','resp_disease', 'viralresp', 'viralgi')
claims_agg<-claims_agg[,(numCols):= lapply(.SD,FUN=function(x) {as.integer(round(as.numeric(x),0))}), .SDcols=numCols]
rm(numCols)

####MERGE IN THE DENOMINATOR####
denom1<-read_sas(paste0(f4,'zip_date_freq_subhosp_',yy,'.sas7bdat'))
#denom1<-read_sas(paste0(f3,'test','.sas7bdat'))
denom1<-setDT(denom1)
denom1[,ZIP_CD:=as.integer(ZIP_CD)]
#Use "all=T" because Robert included all the frequencies for all the dates in all the ZIP codes, so this fills in zeros where no admissions occurred in that ZIP code on that date.
#create a year variable
claims_agg[,year:=as.POSIXlt(CLM_FROM_DT)$year+1900]
#drop if more than 1 year prior (in 2011, this is only 4 rows)
claims_agg<-claims_agg[year>=yy-1,]
dat1<-merge(claims_agg[year==yy,],denom1,by.x=c('ZIP_CD','CLM_FROM_DT'),by.y=c('ZIP_CD','NEWDATE'),all=T)

#for individuals whose claim-from-date is prior to year yy, also import the prior year
dat2<-data.table()
if (yy>2006) {
  v2<-yy-1
  denom1<-read_sas(paste0(f4,'zip_date_freq_subhosp_',v2,'.sas7bdat'))
  denom1<-setDT(denom1)
  denom1[,ZIP_CD:=as.integer(ZIP_CD)]
  dat2<-merge(claims_agg[year==v2,],denom1,by.x=c('ZIP_CD','CLM_FROM_DT'),by.y=c('ZIP_CD','NEWDATE'),all.x=T)
}
claims_zeros<-rbind(dat1,dat2,fill=T)
rm(v2)

#replace NAs for admissions, created with merge with denominator, with zeros, since the NAs now mean no admissions occurred in that ZIP_CD on that date
v1<-c('gi_disease','resp_disease','viralresp','viralgi')
claims_zeros[,v1]<-claims_zeros[,lapply(.SD,FUN=function(x) {x[is.na(x)]<-0;x}),.SDcols=v1]

#drop rows that have a denominator of missing or NA, because we can't have a denominator of zero or missing
claims_zeros<-claims_zeros[!(is.na(NEW_COUNT) | NEW_COUNT==0),]

####MERGE IN THE ZCTA####
claims_zeros<-merge(claims_zeros,zip_to_zcta_mi_oh_pa[,c('ZIP_CD','ZCTA')],by='ZIP_CD',all.x=T)

#identify ZIP codes with missing ZCTAs and drop
v1<-unique(claims_zeros$ZIP_CD[is.na(claims_master$ZCTA)])
#in 2006, ZIP codes 17354, 17583, and 54540 did not have corresponding ZCTAs.
v3<-nrow(claims_red)-nrow(claims_red[ZIP_CD %in% as.character(v1),]) #0 in 2006
claims_dim[year==yy,zip.to.zcta.present:=v3]
#claims_zeros[ZIP_CD %in% v1,mean(COUNT,na.rm=T)] #these represent an average of 4 persons in 2006
claims_zeros<-claims_zeros[!(ZIP_CD %in% v1),]
rm(v1,v3)

####JOIN TO MASTER FILE####
claims_master<-rbindlist(l=list(claims_master,claims_zeros),fill=F,idcol=NULL)
rm(claims_agg,nodup_claims1, claims_zeros, denom1)
save.image(paste0(f1,'GI_Resp_Anal.RData'))
}

####RE-AGGREGATE AND ADD ZEROS####
#Because some of the years had claims-start-dates that preceded that year, these were added on at a different time than their corresponding year's data set. Therefore, the data must be reaggregated.
dat1<-claims_master[,lapply(.SD,sum,na.rm=T),by=c('ZCTA','CLM_FROM_DT','ZIP_CD','NEW_COUNT'),.SDcols=c('gi_disease','resp_disease','viralresp','viralgi')]
#38792 rows had been picked up as admissions starting in a prior year and ending in the year of the corresponding data set.
claims_master<-dat1
rm(dat1,yy)

####IMPORTING PRECIP & MERGING BY ZCTA and DATE####
ppt_master<-data.table()
v2<-claims_master[,unique(ZCTA)]
for (yy in 2006:2013) {
  ppt1<-read_sas(paste0(f3,"zcta13_miohpa_prism_ppt_",yy,".sas7bdat"))
  ppt1<-setDT(ppt1)
  names(ppt1)<-toupper(names(ppt1))
  
  #precipitation is in millimeters, so to save space, we sacrifice precision and convert to whole millimeter values. ZCTA is also converted here to an integer.
  numCols<-c('ZCTA',grep('PRISM',names(ppt1),value=T))
  ppt1<-ppt1[,(numCols):= lapply(.SD,FUN=function(x) {as.integer(round(as.numeric(x),0))}), .SDcols=numCols]
  #checking important variables 
  #summary(ppt1$DATE) #for 2009, 2009-01-01 to 2009-12-31 
  
  cat(yy,'no ppt missing:',all(!is.na(ppt1[,grep('PRISM',names(ppt1))]))) #TRUE, so none missing
  
  #subset to lag days 0-20
  ppt1<-ppt1[,c('ZCTA','DATE','PRISM_PPT',paste0('PRISM_PPT_L',1:20))]
  
  #bind to master
  ppt_master<-rbind(ppt_master,ppt1[ZCTA %in% v2,c('ZCTA','DATE','PRISM_PPT',paste0('PRISM_PPT_L',1:20))])
  rm(ppt1,numCols)
}

#join with claims_master
claims_master<-merge(claims_master,ppt_master,by.x=c('ZCTA','CLM_FROM_DT'),by.y=c('ZCTA','DATE'),all.x=T)
rm(ppt_master)

####IMPORTING TEMPERATURE DATA AND MERGING WITH CLAIMS ####
tmax_master<-data.table()
for (yy in 2006:2013) {
  tmax1<-read_sas(paste0(f3,"zcta13_miohpa_prism_tmax_",yy,".sas7bdat"))
  tmax1<-setDT(tmax1)
  names(tmax1)<-toupper(names(tmax1)) #change all names to uppercase
  
  #precipitation is in millimeters, so to save space, we sacrifice precision and convert to whole millimeter values. ZCTA is also converted here to an integer.
  numCols<-c('ZCTA',grep('PRISM',names(tmax1),value=T)) #identify numeric columns and convert to integer
  tmax1<-tmax1[,(numCols):= lapply(.SD,FUN=function(x) {as.integer(round(as.numeric(x),0))}), .SDcols=numCols]
  
  cat(yy,'no tmax missing:',all(!is.na(tmax1[,grep('PRISM',names(tmax1))]))) #TRUE, so none missing
  
  #subset to lag days 0-20
  tmax1<-tmax1[,c('ZCTA','DATE','PRISM_TMAX',paste0('PRISM_TMAX_L',1:20))]
  
  #bind to tmax_master
  tmax_master<-rbind(tmax_master,tmax1[ZCTA %in% v2,c('ZCTA','DATE','PRISM_TMAX',paste0('PRISM_TMAX_L',1:20))])
  rm(tmax1,numCols)
}

#join to claims_master
claims_master<-merge(claims_master,tmax_master,by.x=c('ZCTA','CLM_FROM_DT'),by.y=c('ZCTA','DATE'),all.x=T)
rm(tmax_master,v2)

####IDENTIFY AND DROP ZCTAs MISSING PRISM####
v1<-claims_master[is.na(PRISM_PPT) | is.na(PRISM_TMAX),unique(ZIP_CD)]
claims_master<-claims_master[!(ZIP_CD %in% v1),]
claims_dim[year==yy,prism.present:=nrow(claims_red[!(ZIP_CD %in% as.character(v1))])]
rm(v1, claims_red)
#claims_master has 6406239 rows

####SAVE####
rm(f2,f3,f4,claimsCols,icdCols,viralgi_cd,viralresp_cd)
#save.image(paste0(f1,'GI_Resp_Anal.RData'))

####EXAMINING DATA####
claims_master[,summary(.SD)]
#no missings. Max precip is 227 mm, which is 57 inches, which is high. Max Tmax is 41 C, which is reasonable. Max daily gi_disease, resp_disease, viralresp, and viralgi, are 16, 20, 6, and 4, which are reasonable. Denominator ranges from 1 to 8282.
claims_master[,hist(PRISM_PPT)]
v1<-claims_master[,log(PRISM_PPT)]
hist(v1)
boxplot(v1)
# histogram of the log of PRISM_PPT looks reasonable and without outliers, although the upper ranges are extremely high levels of precipitation. The maximum recorded precipitation in Pennsylvania was only 24 inches (https://climate.met.psu.edu/data/state/staterecords.php). The maximum recorded precipitation in Ohio was 47 inches (http://coolweather.net/staterainfall/ohio.htm). The maximum recorded precipitation in Michigan was 13 inches (https://www.fox17online.com/2019/09/24/michigan-sets-new-state-rainfall-record/). This higher precipitation may be an artifact of the modeling, and because these are also not outliers, we will retain them.
rm(v1)
test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(gi_disease)]),breaks=seq(0,15))
test1<-claims_master[,hist(gi_disease)]
lines(seq(0.5,15)-0.5,test$counts,col='red')
#GI counts (gray bars) follow a Poisson distribution (red line) except that the counts less than the mean are slightly lower and the counts at the mean are slightly higher.

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(resp_disease)]),breaks=seq(0,20))
test1<-claims_master[,hist(resp_disease)]
lines(seq(0.5,20)-0.5,test$counts,col='red')
#Resp counts are also slightly overdispersed.
#Neither is a true Poisson distribution, although neither needs zero-inflation methods either as the number of zeros conforms to a standard Poisson distribution.
rm(test,test1)

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(viralresp)]),breaks=seq(0,6,by=0.5))
test1<-claims_master[,hist(viralresp,breaks=seq(0,6,by=0.5))]
lines(seq(0.5,6,by=0.5)-0.25,test$counts,col='red')
#This is a Poisson distribution.

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(viralgi)]),breaks=seq(0,4,by=0.5))
test1<-claims_master[,hist(viralgi,breaks=seq(0,4,by=0.5))]
lines(seq(0.5,4,by=0.5)-0.25,test$counts,col='red')
#This is also a Poisson distribution.


