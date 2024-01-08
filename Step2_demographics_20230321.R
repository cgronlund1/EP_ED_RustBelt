#################### Step 2 Demographic table code  ####################
####################             Madeline          #####################
## Updated 2022-11-27 by Carina                          ##
#################################################################################
################################################################################

####LIBRARIES####
source('T:\\Analysis\\Madeline\\MS_Precip_Libraries.R')

####LOAD FILES####
#load(paste0(f1,'temp.RData')) #this is added onto in step 3 program, so be careful in modifying it
#load(paste0(f1,'GI_Resp_Anal.RData'))

#Note: in this temp.RData workspace, the variable names are PRISM_PPT_xx, but these are actually ln(PRISM_PPT).

####ADDITIONAL DATA MANAGEMENT####
#load county-to-zip crosswalk
cnty_to_zip<-read.csv(paste0(f5,'zcta2cntys2010geocorr.csv'),stringsAsFactors=F,skip=1)
names(cnty_to_zip)[1]<-'ZCTA'

#keep the county assignment that has the highest population
cnty_to_zip<-cnty_to_zip[order(cnty_to_zip$ZCTA,cnty_to_zip$zcta5.to.county.alloc.factor,decreasing=T),]
cnty_to_zip<-cnty_to_zip[!duplicated(cnty_to_zip$ZCTA),]

precipCols<-c('PRISM_PPT',paste0('PRISM_PPT_L', 1:20))
tempCols<-c('PRISM_TMAX',paste0('PRISM_TMAX_L', 1:20))

#to reduce size, identify ZIP codes where the offset is ever less than 100 and drop these
v1<-claims_master[NEW_COUNT<100,unique(ZIP_CD)]
claims_master<-claims_master[!(claims_master$ZIP_CD %in% v1),]
rm(v1)
#claims_master goes from 6,406,239 to 4,534,579 rows

#merge in county
claims_master<-merge(claims_master,cnty_to_zip[,c('ZCTA','county')],by='ZCTA',all.x=T)
#claims_master[,table(county,useNA='a')] no NAs

#To make the data set small enough for the models to run, keep only the 10 counties in each state with the highest number of homes built before 1940. From American Community Survey 2016 5-year estimates, the 10 counties in each of MI, OH, and PA with the highest number of homes built before 1940 are: Michigan counties of Genesee (26049), Ingham (26065), Kent (26081), Oakland (26125), Washtenaw (26161), Macomb (26099), Kalamazoo (26077), Saginaw (26145), Calhoun (26025) and Wayne (26163), the Ohio counties of Cuyahoga (39035), Lucas (39095), Hamilton (39061), Franklin (39049), Montgomery (39113), Stark (39151), Mahoning (39099), Lorain (39093), Butler (39017) and Summit (39153), and the Pennsylvania counties of Allegheny (42003), Delaware (42045), Luzerne (42079), Montgomery (42091), Berks (42011), Lancaster(42071), Lackawanna (42069), Westmoreland (42129), York (42133) and Philadelphia (42101).
v1<-c(26049,26065,26081,26125,26161,26099,26077,26145,26025,26163,39035,39095,39061,39049,39113,39151,39099,39093,39017,39153,42003,42045,42079,42091,42011,42071,42069,42129,42133,42101) #added 26049 on 210428 because it was left off on accident in 210408 analyses
claims_master<-claims_master[county %in% v1,]
rm(v1)
#2810599 rows on 21-04-08 and 2874883 rows on 210428 when including 26049
claims_master[,sum(gi_disease,na.rm=T)] #2375844
claims_master[,sum(resp_disease,na.rm=T)] #2715456

claims_master[,table(county,useNA='a')][order(claims_master[,table(county,useNA='a')])]

#Log transforming precipitation-- converting all 0 and neg numbers to 1 mm and adding 1 mm to everything else
lnppt1<-claims_master[,lapply(.SD,function(x) ifelse(x <= 0, log(1), log(x + 1))),.SDcols=precipCols]
#multiply by 100 and convert to matrix of integers to save space
lnppt1<-as.matrix(lnppt1[,lapply(.SD,FUN=function(x) {as.integer(round(x*100,0))})])
#to back transform to mm, first divide by 100, then take exp, then subtract 1

v1<-grep('PRISM_PPT',names(claims_master),invert=T)
claims_master<-claims_master[,..v1]
claims_master<-cbind(claims_master,lnppt1)
rm(lnppt1)

##Identify ZIP_CDs that have the same PRISM_PPT values (actually ln(PPT)*100+1) values when rounded to the nearest tens place. This will help reduce spatial autocorrelation, as these are often adjacent.
#First, find the unique PRISM series.
v1<-grep('PPT',names(claims_master),value=T)
test<-claims_master[,.SD,.SDcols=c('CLM_FROM_DT','county','ZIP_CD',v1)]
test[,(v1):=lapply(.SD,FUN=function(x) {as.integer(round(x/10,0))}),.SDcols=v1]
test1<-dcast(test[,],as.formula(c('county+ZIP_CD~CLM_FROM_DT')),value.var=as.list(v1))
test1[,PRISMseries:=apply(test1[,3:61364],1,FUN=paste,collapse='_')]
uniqueN(test1[,PRISMseries]) #824 unique series
test1<-test1[,c('county','ZIP_CD','PRISMseries')]
#Label the unique series as the first county and first ZIP code with that series value with a 9 in front of it and merge it to the ZIP codes to indicate which ZIP code goes with which series, or super ZIP.
test1<-test1[order(county,ZIP_CD),]
test2<-unique(test1,by='PRISMseries')
test2[,super.zip:=paste0('s',county,ZIP_CD)]
test1<-test2[,c('PRISMseries','super.zip')][test1,on='PRISMseries']
super.zips<-test1[,c('ZIP_CD','county','super.zip')]
super.zips[,N:=.N,by='super.zip']
super.zips[,table(N)]
#  1   2   3   4   8 
#711 158  75  32   8
#So there are 158 instances of 2 ZIP codes sharing a PRISM series, 75 instances of 3 ZIP codes sharing a PRISM series, etc.
rm(test,test1,test2,v1)

#Merge the super ZIP names into the master claims file by county and ZIP code.
claims_master<-super.zips[,c('county','ZIP_CD','super.zip')][claims_master,on=c('county','ZIP_CD')]

#plot each of the ZIP_CD sets of unique series
load(paste0(f1,'study_zctas_shp.RData'))
for (i in unique(super.zips[N>1,super.zip])) {
  cat(i)
  v2<-super.zips[super.zip==i,as.character(ZIP_CD)]
  v1<-unique(substr(v2,1,2))
  plot(study.zctas.shp[substr(study.zctas.shp$ZCTA5CE10,1,2) %in% v1,])
  plot(study.zctas.shp[study.zctas.shp$ZCTA5CE10 %in% v2,],add=T,col='red')
  Sys.sleep(1)
}
#In looking at these 273 plots, these were all instances of ZIP codes being adjacent to each other and likely sharing a PRISM grid cell. If we were to further aggregate precipitation and temperature, e.g., by 10s of degrees C and cm instead of mm, or if we didn't use all 21 days, then we might see instances of non-adjacent ZIP codes a sharing PRISM series or of ZIP codes that were not in the same grid cell but were still close sharing a PRISM series.

#for each new super ZCTA, find its centroid
#read in the ZCTA centroids
zcta.centroids<-fread('T:\\Data Management\\zipzctaxwalks\\zcta_popwt_centroids.csv')
zcta.centroids<-zcta.centroids[ZCTA %in% claims_master[,ZCTA],]
#merge ZIP code to ZCTA to centroid and take mean by super ZIP
dt1<-setDT(zip_to_zcta_mi_oh_pa)[,.(ZIP_CD,ZCTA)][super.zips,on='ZIP_CD']
dt2<-zcta.centroids[,.(ZCTA,XCoord,YCoord,POINT_X,POINT_Y)][dt1,on='ZCTA']
v1<-c('XCoord','YCoord','POINT_X','POINT_Y')
dt2[,paste0(v1,'.mean'):=lapply(.SD,mean),.SDcols=v1,by='super.zip']
super.zips<-dt2
rm(dt1,dt2,v1)

#Further reduce claims_master by super ZIPs. This helps with spatial autocorrelation and slightly reduces computation time. This takes a few minutes to run.
v1<-grep('PRISM',names(claims_master),value=T)
dt1<-claims_master[,lapply(.SD,FUN=function(x) {as.integer(sum(x))}),.SDcols=c('NEW_COUNT','gi_disease','resp_disease','viralresp','viralgi'),by=c('CLM_FROM_DT','super.zip')]
dt2<-claims_master[,lapply(.SD,FUN=function(x) {as.integer(round(mean(x),0))}),.SDcols=v1,by=c('CLM_FROM_DT','super.zip')]
all(dt1[,1:2]==dt2[,1:2]) #TRUE, so can do cbind instead of merge
claims_master_agg<-cbind(dt1,dt2[,3:44])
rm(dt1,dt2)

####create spatial variables####
#add county variable back to claims_master_agg. When there is more than one county for a super ZIP, we assign the county in which most of the super ZIP's population resides.
#first, merge the ZCTA population info to the ZIP codes
dt1<-setDT(cnty_to_zip)
dt2<-zip_to_zcta_mi_oh_pa[,.(ZIP_CD,ZCTA)][dt1,on='ZCTA']
dt2<-dt2[!is.na(ZIP_CD),]
#next, merge in the super zips
dt2<-super.zips[,.(ZIP_CD,super.zip)][dt2,on='ZIP_CD']
dt3<-dt2[!is.na(super.zip),sum(Total.Pop..2010.census),by=.(super.zip,county)]
dt3[,table(super.zip)][dt3[,table(super.zip)>1]]
#Super zips s2609948080,s2609948092,s2612548167,s4204519023,s4204519041,s4204519050,s4209119012, and s4209119038 all lie across two counties.
#next, select the super.zip-county combinations for each super.zip with the largest population and use that county
dt4<-dt3[,.SD[which.max(V1)],by=super.zip]
dt5<-dt4[,.(super.zip,county)]
#additionally, merge the mean X and Y coordinates, which serve as the centroids for the super.zips, into the superzip_to_1cnty file which also has only unique values for super.zips.
dt6<-unique(super.zips[,.(super.zip,XCoord.mean,YCoord.mean,POINT_X.mean,POINT_Y.mean)])
superzip_to_1cnty<-dt6[dt5,on='super.zip']

####save because it takes a few minutes to generate these data####
save(list=c('claims_dim','claims_master_agg','cnty_to_zip','colNamesSAS','study.zctas.shp','super.zips','zcta.centroids','zip_to_zcta_mi_oh_pa','precipCols','tempCols','superzip_to_1cnty'),file=paste0(f1,'temp.RData'))

####create and save summaries of final data####
v1<-c(26065,26081,26125,26161,26099,26077,26145,26025,26163,39035,39095,39061,39049,39113,39151,39099,39093,39017,39153,42003,42045,42079,42091,42011,42071,42069,42129,42133,42101)
v2<-claims_master[,unique(ZIP_CD)]

claims_master[,length(unique(super.zip))] #824 super ZIPs
do.call(rbind,lapply(claims_master[,.(PRISM_PPT/10, PRISM_TMAX)],FUN=quantile,c(0,0.75,0.9,0.95,0.995,1)))
#            0% 75% 90%  95% 99.5% 100%
#V1           0  11  24 28.9  38.1 54.2
#PRISM_TMAX -17  26  30 31.0  35.0 41.0

####Plot the precipitation time series####
claims_master<-cbind(claims_master,lnppt1)
v1<-claims_master[,quantile(exp(PRISM_PPT/100)-1,0.75),by=as.POSIXlt(CLM_FROM_DT)$yday]
v2<-claims_master[,exp(max(PRISM_PPT/100))-1,by=as.POSIXlt(CLM_FROM_DT)$yday]
v3<-claims_master[,exp(min(PRISM_PPT/100))-1,by=as.POSIXlt(CLM_FROM_DT)$yday]
v4<-claims_master[,exp(quantile(PRISM_PPT/100,0.95))-1,by=as.POSIXlt(CLM_FROM_DT)$yday]
plot(v1,type='l',ylim=c(0,max(v2[,2])),xlab=c('Day of Year'),ylab=c('Precipitation (mm)'))
lines(v2,col='red')
lines(v3,col='red')
lines(v4,col='green')
legend('topleft',lty=1,col=c('black','green','red'),legend=c('75th percentile','95th percentile','Minimum and Maximum'))

####Plot the precipitation histograms####
claims_master[,hist(exp(PRISM_PPT/100)-1,xlab='Precipitation (mm)',main='')]

claims_master[,hist(PRISM_PPT/100,xlab='ln(Precipitation)',main='')]

####Histograms of outcomes####
#These are the daily ED visits for a specific cause in a ZIP code across 8 years
par(mfrow=c(1,2))
test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(gi_disease)]),breaks=seq(0,15),plot=F)
test1<-claims_master[,hist(gi_disease, main='',xlab='Gastrointestinal Disease Visits')]
lines(seq(0.5,15)-0.5,test$counts,col='red')
#GI counts (gray bars) follow a Poisson distribution (red line) except that the counts less than the mean are slightly lower and the counts at the mean are slightly higher.

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(resp_disease)]),breaks=seq(0,20),plot=F)
test1<-claims_master[,hist(resp_disease,main='',xlab='Respiratory Disease Visits')]
lines(seq(0.5,20)-0.5,test$counts,col='red')
#Resp counts are also slightly overdispersed.
#Neither is a true Poisson distribution, although neither needs zero-inflation methods either as the number of zeros conforms to a standard Poisson distribution.
rm(test,test1)

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(viralresp)]),breaks=seq(0,6,by=0.5),plot=F)
test1<-claims_master[,hist(viralresp,breaks=seq(0,6,by=0.5))]
lines(seq(0.5,6,by=0.5)-0.25,test$counts,col='red')
#This is a Poisson distribution.

test<-hist(rpois(n=nrow(claims_master),lambda=claims_master[,mean(viralgi)]),breaks=seq(0,4,by=0.5),plot=F)
test1<-claims_master[,hist(viralgi,breaks=seq(0,4,by=0.5))]
lines(seq(0.5,4,by=0.5)-0.25,test$counts,col='red')
#This is also a Poisson distribution.

####ED VISIT TIME SERIES####
#These are the daily ED visits for GI disease in a ZIP code across 8 years
par(mfrow=c(1,1))
v1<-claims_master[,mean(gi_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v2<-claims_master[,max(gi_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v3<-claims_master[,min(gi_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v4<-claims_master[,quantile(gi_disease/NEW_COUNT,0.95),by=as.POSIXlt(CLM_FROM_DT)$yday]
plot(v1,type='l',ylim=c(0,max(v2[,2])*1.1),xlab=c('Day of Year'),ylab=c('GI ED Visits'))
lines(v2,col='red')
lines(v3,col='red')
lines(v4,col='blue')
legend('topright',lty=1,col=c('black','blue','red'),legend=c('Mean','95th percentile','Minimum and Maximum'),cex=.8,ncol=3)

#These are the daily ED visits for Respiratory disease in a ZIP code across 8 years
par(mfrow=c(1,1))
v1<-claims_master[,mean(resp_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v2<-claims_master[,max(resp_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v3<-claims_master[,min(resp_disease),by=as.POSIXlt(CLM_FROM_DT)$yday]
v4<-claims_master[,quantile(resp_disease,0.95),by=as.POSIXlt(CLM_FROM_DT)$yday]
plot(v1,type='l',ylim=c(0,max(v2[,2])+10),xlab=c('Day of Year'),ylab=c('Respiratory ED Visits'))
lines(v2,col='red')
lines(v3,col='red')
lines(v4,col='blue')
legend('topright',lty=1,col=c('black','blue','red'),legend=c('Mean','95th percentile','Minimum and Maximum'),cex=0.8,ncol=3)

####MEAN DAILY ED VISIT RATE BY YEAR####
par(mfrow=c(1,2))
claims_master[,year:=1900+year]
claims_master[,mean(gi_disease/NEW_COUNT)*100000,by=year]
#year       V1
#1: 2006 54.96794
#2: 2007 56.09813
#3: 2008 57.75552
#4: 2009 58.04373
#5: 2010 59.09233
#6: 2011 61.46405
#7: 2012 62.38638
#8: 2013 63.12102

claims_master[,mean(resp_disease/NEW_COUNT)*100000,by=year]
#  year       V1
#1: 2006 67.24003
#2: 2007 66.13340
#3: 2008 67.76065
#4: 2009 66.43726
#5: 2010 66.73093
#6: 2011 69.14931
#7: 2012 68.79224
#8: 2013 70.55549

#ED visit rate distributions
t(claims_master_agg[,lapply(.SD,FUN=function(x) round(quantile(x/NEW_COUNT*10000,p=c(0,0.25,0.5,0.75,0.9,0.95,0.99,0.995,0.999,1)),0)),.SDcols=c('gi_disease','resp_disease')])
#             [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#gi_disease      0    0    3    9   15   21   41   53   78   268
#resp_disease    0    0    4   10   17   23   45   57   82   263

####EXTREME PRECIPITATION EVENTS####
#Extreme precipitation events 
#defined as days with precipitation in the top 1 percent of all days with precipitation.
claims_master[,quantile(PRISM_PPT,c(0.5,0.75,0.9,0.95,.99,.995,.999))]
# 50%   75%   90%   95%   99% 99.5% 99.9% 
#   0   110   240   289   358   383   429

claims_master[,quantile(exp(PRISM_PPT/100)-1,c(0.5,0.75,0.9,0.95,.99,.995,.999))]
#50%       75%       90%       95%       99%     99.5%     99.9% 
#0.000000  2.004166 10.023176 16.993310 34.873541 45.062538 71.966468 

nrow(claims_master[PRISM_PPT>=358,]) #29,989 ZCTA-days from 2006-2013
claims_master[PRISM_PPT>=358,sum(gi_disease)]/claims_master[,sum(gi_disease)]
#risk of 0.00992 EP events per person, which corresponds to about 1%

#Temps in cold season
claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% c(10,11,0,1),quantile(PRISM_TMAX,c(0,0.25,0.5,0.75,1))]

####DEMOGRAPHICS USING CLAIMS LEVEL DATA#####
#See Robert's summaries from SAS and also calculate from the scrambled data

####Maps####
library(tigris)
library(sp)
library(rgdal)
library(RColorBrewer)

studarea<-counties(state=c('MI','OH','PA'),cb=T)
plot(studarea['GEOID'],col='white')

countyBens<-claims_master_agg[,sum(NEW_COUNT)/(365*8),by=substr(super.zip,2,6)]
#countyPrecip<-claims_master_agg[,mean((exp(PRISM_PPT/100)-1)),by=substr(super.zip,2,6)] #mean
countyPrecip<-claims_master_agg[,exp(mean(PRISM_PPT/100))-1,by=substr(super.zip,2,6)] #geometric mean
#countyPrecip<-claims_master_agg[,mean(exp(PRISM_PPT/100)-1>25.4),by=substr(super.zip,2,6)] #number of days > 1 in
setnames(countyBens,c('substr','V1'),c('GEOID','meanDailyBens'))
setnames(countyPrecip,c('substr','V1'),c('GEOID','meanDailyPrecip'))
dat1<-merge(studarea,countyBens,by='GEOID',all.x=T)
dat1<-merge(dat1,countyPrecip,by='GEOID',all.x=T)
brewer.pal.info

v1<-colorRampPalette(colors=c('orange','darkred'))(8)
v4<-colorRampPalette(colors=c('orange','darkred'))(6)
quantile(dat1$meanDailyBens,p=seq(0,1,length.out=9),na.rm=T)
#       0%     12.5%       25%     37.5%       50%     62.5%       75%     87.5%      100% 
#14502.95  20084.02  22952.39  30425.52  38444.00  41591.29  57498.41  70587.01 147503.42 
dat1$cols<-as.numeric(cut(dat1$meanDailyBens,breaks=c(14500,20000,30000,40000,50000,70000,148000)))
quantile(dat1$meanDailyPrecip,p=seq(0,1,length.out=3),na.rm=T)
#       0%       50%      100% 
#0.7617002 0.9071049 1.1599448
dat1$cols.precip<-round(dat1$meanDailyPrecip-0.7617,1)
v2<-colorRampPalette(colors = c('lightblue','navyblue'))(13)
v3<-colorRampPalette(colors=c('lightblue','navyblue'))(13)
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))
#layout(mat=matrix(c(1,2,3),1,3))
#layout.show()
plot(dat1['GEOID'],col=v1[dat1$cols],main='(A)',adj=0)
legend('topright',legend=c('15k-20k','20k-30k','30k-40k','40k-50k','50k-70k','70k-148k'),col=v4,fill=v4,cex=1,title='Beneficiaries')
#text(x=1,y=1,'(A)',adj=c(1,1))
plot(dat1['GEOID'],col=v2[dat1$cols.precip*10+1],main='(B)',adj=0)
legend('topright',legend=c(rep(NA,2),'0.8 mm',rep(NA,5),'0.9 mm',rep(NA,5),'1.2 mm'),fill=c('white','white',v2),y.intersp=0.5,border=NA,cex=1,title='Precipitation',title.adj=0.5,xpd=T)
text(x=-85,y=48.5,'(B)',adj=c(1,1))

####Precipitation Waves####
#claims_master_agg[,quantile((exp(PRISM_PPT/100)-1)/10,c(0.75,0.9,0.95,0.99))]
#75%       90%       95%       99% 
#0.2004166 1.0023176 1.6993310 3.4873541
#In cm, the 90th percentile is 1 cm, the 95th is 1.7, and the 99th is 3.5.
#How often are there 2 2-cm days in a row?
claims_master_agg[,p2cm0_1:=(exp(PRISM_PPT/100)-1)/10>=2 & (exp(PRISM_PPT_L1/100)-1)/10>=2][,p2cm0_1:=as.integer((p2cm0_1-mean(p2cm0_1))*10000)]
claims_master_agg[,table(p2cm0_1)/length(p2cm0_1)*100]
#only 0.34% of the time

#How often are there 3 2-cm days in a row?
claims_master_agg[,p2cm0_2:=(exp(PRISM_PPT/100)-1)/10>=2 & (exp(PRISM_PPT_L1/100)-1)/10>=2 & (exp(PRISM_PPT_L2/100)-1)/10>=2][,p2cm0_2:=as.integer((p2cm0_2-mean(p2cm0_2))*10000)]
claims_master_agg[,table(p2cm0_2)/length(p2cm0_2)*100]
#0.035% of the time, or only once every 8 years, which is too infrequent for this study's time period of 8 years
claims_master_agg[,p2cm0_2:=NULL]

#How often are there 2 1-cm days in a row?
claims_master_agg[,p1cm0_1:=as.integer((exp(PRISM_PPT/100)-1)/10>=1 & (exp(PRISM_PPT_L1/100)-1)/10>=1)][,p1cm0_1:=as.integer((p1cm0_1-mean(p1cm0_1))*10000)]
claims_master_agg[,table(p1cm0_1)/length(p1cm0_1)*100]
#1.8% of the time

#How often are there 3 1-cm days in a row?
claims_master_agg[,p1cm0_2:=as.integer((exp(PRISM_PPT/100)-1)/10>=1 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10>=1)][,p1cm0_2:=as.integer((p1cm0_2-mean(p1cm0_2))*10000)]
claims_master_agg[,table(p1cm0_2)/length(p1cm0_2)*100]
#0.27% of the time, or once a year

#How often are there 4 1-cm days in a row?
claims_master_agg[,p1cm0_3:=as.integer((exp(PRISM_PPT/100)-1)/10>=1 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10>=1 & (exp(PRISM_PPT_L3/100)-1)/10>=1)][,p1cm0_3:=as.integer((p1cm0_3-mean(p1cm0_3))*10000)]
claims_master_agg[,table(p1cm0_3)/length(p1cm0_3)*100]
#0.044% of the time, or once every 6 years, which is too infrequent for this study's time period
claims_master_agg[,p1cm0_3:=NULL]

#How often are there 5 0-cm days in a row?
claims_master_agg[,p0cm0_4:=as.integer((exp(PRISM_PPT/100)-1)/10==0 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10==0 & (exp(PRISM_PPT_L3/100)-1)/10==0 & (exp(PRISM_PPT_L4/100)-1)/10==0)][,p0cm0_4:=as.integer((p0cm0_4-mean(p0cm0_4))*10000)]
claims_master_agg[,table(p0cm0_4)/length(p0cm0_4)*100]
#1.12% of the time

#How often are there 6 0-cm days in a row?
claims_master_agg[,p0cm0_5:=as.integer((exp(PRISM_PPT/100)-1)/10==0 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10==0 & (exp(PRISM_PPT_L3/100)-1)/10==0 & (exp(PRISM_PPT_L4/100)-1)/10==0 & (exp(PRISM_PPT_L5/100)-1)/10==0)][,p0cm0_5:=as.integer((p0cm0_5-mean(p0cm0_5))*10000)]
claims_master_agg[,table(p0cm0_5)/length(p0cm0_5)*100]
#0.8% of the time

##How often are there 7 0-cm days in a row?
claims_master_agg[,p0cm0_6:=as.integer((exp(PRISM_PPT/100)-1)/10==0 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10==0 & (exp(PRISM_PPT_L3/100)-1)/10==0 & (exp(PRISM_PPT_L4/100)-1)/10==0 & (exp(PRISM_PPT_L5/100)-1)/10==0 & (exp(PRISM_PPT_L6/100)-1)/10==0)][,p0cm0_6:=as.integer((p0cm0_6-mean(p0cm0_6))*10000)]
claims_master_agg[,table(p0cm0_6)/length(p0cm0_6)*100]
#0.57% of the time

##How often are there 8 0-cm days in a row?
claims_master_agg[,p0cm0_7:=as.integer((exp(PRISM_PPT/100)-1)/10==0 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10==0 & (exp(PRISM_PPT_L3/100)-1)/10==0 & (exp(PRISM_PPT_L4/100)-1)/10==0 & (exp(PRISM_PPT_L5/100)-1)/10==0 & (exp(PRISM_PPT_L6/100)-1)/10==0 & (exp(PRISM_PPT_L7/100)-1)/10==0)][,p0cm0_7:=as.integer((p0cm0_7-mean(p0cm0_7))*10000)]
claims_master_agg[,table(p0cm0_7)/length(p0cm0_7)*100]
#0.4% of the time

##How often are there 9 0-cm days in a row?
claims_master_agg[,p0cm0_8:=as.integer((exp(PRISM_PPT/100)-1)/10==0 & (exp(PRISM_PPT_L1/100)-1)/10>=1 & (exp(PRISM_PPT_L2/100)-1)/10==0 & (exp(PRISM_PPT_L3/100)-1)/10==0 & (exp(PRISM_PPT_L4/100)-1)/10==0 & (exp(PRISM_PPT_L5/100)-1)/10==0 & (exp(PRISM_PPT_L6/100)-1)/10==0 & (exp(PRISM_PPT_L7/100)-1)/10==0 & (exp(PRISM_PPT_L8/100)-1)/10==0)][,p0cm0_8:=as.integer((p0cm0_8-mean(p0cm0_8))*10000)]
claims_master_agg[,table(p0cm0_8)/length(p0cm0_8)*100]
#0.28% of the time

#How often are there two 1/2 inch days in a row?
claims_master_agg[,p05in0_1:=as.integer((exp(PRISM_PPT/100)-1)/10/2.54>=0.5 & (exp(PRISM_PPT_L1/100)-1)/10/2.54>=0.5)][,p05in0_1:=as.integer((p05in0_1-mean(p05in0_1))*10000)]
claims_master_agg[,table(p05in0_1)/length(p05in0_1)*100]
#1.1% of the time

#How often are there three 1/2 inch days in a row?
claims_master_agg[,p05in0_2:=as.integer((exp(PRISM_PPT/100)-1)/10/2.54>=0.5 & (exp(PRISM_PPT_L1/100)-1)/10/2.54>=0.5 & (exp(PRISM_PPT_L2/100)-1)/10/2.54>=0.5)][,p05in0_2:=as.integer((p05in0_2-mean(p05in0_2))*10000)]
claims_master_agg[,table(p05in0_2)/length(p05in0_2)*100]
#0.14% of the time, or 4 times over 8 years
claims_master_agg[,table(p05in0_2,as.POSIXlt(CLM_FROM_DT)$mon)]
#p05in0_2      0      1      2      3      4      5      6      7      8      9     10     11
#-13  204286 186160 204138 197343 204180 197324 203975 203787 197006 203987 197577 204314
#9986     35     36    183    387    141    406    346    534    724    334    153      7
#So only 7 zip-code-days in December.

#How often are there two 1 inch days in a row?
claims_master_agg[,p1in0_1:=as.integer((exp(PRISM_PPT/100)-1)/10/2.54>=1 & (exp(PRISM_PPT_L1/100)-1)/10/2.54>=1)][,p1in0_1:=as.integer((p1in0_1-mean(p1in0_1))*10000)]
claims_master_agg[,table(p1in0_1)/length(p1in0_1)*100]
#0.17% of the time, or 5 times over 8 years
claims_master_agg[,table(p1in0_1,as.POSIXlt(CLM_FROM_DT)$mon)]
#p1in0_1      0      1      2      3      4      5      6      7      8      9     10     11
#-17  204318 185999 204207 197377 204225 197023 203902 203799 196864 203753 197568 204232
#9982      3    197    114    353     96    707    419    522    866    568    162     89
#Only 3 ZIP-code-days in January


####END####
