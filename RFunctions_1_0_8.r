########################################################################################################
#R code for diagnose, balance, lengthen, makeplot, omit.history, makehistory.one, makehistory.two, 
#and apply.scope  functions                                                
#
#These functions can be used to implement the diagnostic framework outlined in:
#Jackson JW. Diagnostics for confounding of time-varying and other joint effects. Epidemiology. 2016.
#
#© John W. Jackson 2015
#THE STANDARD MIT LICENSE APPLIES:
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software
#and associated documentation files (the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
#and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial
#portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
#LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A #PARTICULAR PURPOSE AND NONINFRINGEMENT.
#IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION #OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
########################################################################################################

##########################
##LOAD REQUIRED PACKAGES##
##########################

library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(Rmpfr)
library(data.table)

###########################
##MAKEHISTORY() FUNCTIONS##
###########################

#Function to create exposure history for a single time varying exposure
makehistory.one <- function (input,times,group=NULL,exposure,name.history="h") {

  if (is.null(input)) {
    stop ("ERROR: 'input' dataframe is missing")
  }

  if (is.null(exposure)) {
    stop ("ERROR: root name for exposure is missing")
  }

  if (is.null(times)) {
    stop ("ERROR: indices for exposure measurement times is missing. Please specify a numeric vector of times")
  }

  list.exposure <- paste(exposure,times,sep="_")
  history.names <- paste(name.history,times,sep="_")

  if (any(!list.exposure %in% names(input))) {
    stop("ERROR: The exposure root name is misspelled, or some exposure measurements are missing from the input dataframe, or incorrect measurement times have been specified")
  }

  matrix.exposure <- matrix(0,nrow=nrow(input),ncol=length(list.exposure))
  matrix.exposure <-input[,list.exposure]
  matrix.history <- matrix(0,nrow=nrow(input),ncol=length(list.exposure))
  colnames(matrix.history) <- history.names
  for (i in 1:length(list.exposure)) {
    if (i==1) {
      matrix.history[,i]<-paste("H",input[,group],sep="")
    } else if (i>1) {
      matrix.history[,i] <- paste(matrix.history[,(i-1)],matrix.exposure[,i-1],sep="")
    }
  }
  data.frame(input,matrix.history)

}

#Function to create joint exposure history for two distinct time-varying exposures
makehistory.two <- function (input,group=NULL,exposure.a,exposure.b,name.history.a="ha",name.history.b="hb",times) {

  list.exposure.a <- paste(exposure.a,times,sep="_")
  list.exposure.b <- paste(exposure.b,times,sep="_")
  history.names.a <- paste(name.history.a,times,sep="_")
  history.names.b <- paste(name.history.b,times,sep="_")

  if (is.null(input)) {
    stop ("ERROR: 'input' dataframe is missing")
  }

  if (is.null(exposure.a)) {
    stop ("ERROR: root name for the first exposure is missing")
  }

  if (is.null(exposure.b)) {
    stop ("ERROR: root name for the second exposure is missing")
  }

  if (is.null(times)) {
    stop ("ERROR: indices for exposure measurement times is missing. Please specify a numeric vector of times")
  }

  if (any(!list.exposure.a %in% names(input))) {
    stop("ERROR: The exposure root name is misspelled, or some exposure measurements are missing from the input dataframe, or incorrect measurement times have been specified")
  }

  if (any(!list.exposure.b %in% names(input))) {
    stop("ERROR: The exposure root name is misspelled, or some exposure measurements are missing from the input dataframe, or incorrect measurement times have been specified")
  }

  matrix.exposure.a <- matrix(0,nrow=nrow(input),ncol=length(list.exposure.a))
  matrix.exposure.b <- matrix(0,nrow=nrow(input),ncol=length(list.exposure.b))

  matrix.exposure.a <- input[,list.exposure.a]
  matrix.exposure.b <- input[,list.exposure.b]

  matrix.history.a <- matrix(0,nrow=nrow(input),ncol=length(list.exposure.a))
  matrix.history.b <- matrix(0,nrow=nrow(input),ncol=length(list.exposure.b))

  joint.exposure <- matrix(0,nrow=nrow(input),ncol=max(length(list.exposure.a),length(list.exposure.b)))

  colnames(matrix.history.a) <-  history.names.a
  colnames(matrix.history.b) <-  history.names.b
  colnames(joint.exposure)    <- paste("q",seq(from=0,to=(max(length(list.exposure.a),length(list.exposure.b))-1)),sep="_")

  for (i in 1:max(length(list.exposure.a),length(list.exposure.b))) {
    if (i==1) {
      matrix.history.a[,i]<-paste("H",input[,group],sep="")
      matrix.history.b[,i]<-as.matrix(paste("H",input[,group],matrix.exposure.a[,i],sep=""))
      joint.exposure[,i] <- as.matrix(paste("H",input[,group],matrix.exposure.a[,i],matrix.exposure.b[,i],sep=""))

    } else if (i>1) {
      matrix.history.a[,i] <- joint.exposure[,(i-1)]
      matrix.history.b[,i] <- paste(joint.exposure[,(i-1)],matrix.exposure.a[,i],sep="")
      joint.exposure[,i] <- as.matrix(paste(joint.exposure[,i-1],matrix.exposure.a[,i],matrix.exposure.b[,i],sep=""))
    }
  }
  data.frame(input,matrix.history.a,matrix.history.b)
}



#######################
##LENGTHEN() FUNCTION##
#######################

#Function to create a "tidy" dataframe where the key observation is the paring of exposure and covariate measurement times
lengthen <- function (input,
                      diagnostic,
                      censoring,
                      id,
                      times.exposure,
					  times.covariate,
                      exposure,
                      temporal.covariate,
                      static.covariate=NULL,
                      history=NULL,
                      weight.exposure=NULL,
                      censor=NULL,
                      weight.censor=NULL,
                      strata=NULL) {

  if(is.null(input)) {
    stop ("ERROR: 'input' is missing or misspecified. Please specify a dataframe in 'wide' format")
  }
  if(is.null(id)) {
    stop ("ERROR: 'id' is missing or misspecified. Please specify a unique identifier for each observation")
  }
  if (!all(!input$id %in% input$id[duplicated(input$id)])) { 
    stop("ERROR: id does not uniquely identify each observation (i.e. each row). Please specify a unique identifier.")
  }
  if(is.null(diagnostic) | !diagnostic %in% c(1,2,3)) {
    stop ("ERROR: 'diagnostic' is missing or misspecified. Please specify as 1, 2 or 3")
  }
  if(is.null(censoring) | !censoring %in% c("no","yes")) {
    stop ("ERROR: 'censoring' is missing. Please specify it as yes or no")
  }
  if(is.null(exposure)) {
    stop ("ERROR: 'exposure' is missing. Please specify the root name for exposure")
  }
  if(is.null(times.exposure) | is.null(times.covariate)) {
    stop ("ERROR: either 'times.exposure' or 'times.covariate' is missing. Please specify an integer  or a numeric vector of times")
  }
  if(is.null(temporal.covariate) & is.null(static.covariate)) {
    stop ("ERROR: both 'temporal.covariate' and 'static.covariate' are missing. Please specify a character vector of root names for temporal.covariates or static.covariates")
  }

  if (censoring=="yes" & is.null(censor)) {
    stop ("ERROR: 'censor' is missing. Please specify a root name for censoring indicators")
  }

  if (diagnostic==1) {
    if (is.null(history)) stop ("ERROR: 'history' is missing. Please specify a root name for exposure history")
  }  else if ((diagnostic==2 | diagnostic==3) & ((is.null(history) | is.null(weight.exposure)) & is.null(strata))) {
    stop ("ERROR: Please specify the root names for exposure history and exposure weight variables, or the root name for strata variables")
  }

  list.exposure  <- paste(exposure,times.exposure,sep="_")

if (!is.null(history))  {
  list.history   <- paste(history,times.exposure,sep="_")
} else {list.history <- NULL
}

if (!is.null(weight.exposure))  {
  list.weight.exposure   <- paste(weight.exposure,times.exposure,sep="_")
} else {list.weight.exposure <- NULL
}

if (!is.null(censor) & diagnostic!=2)  {
  list.censor   <- paste(censor,times.exposure,sep="_")
} else {list.censor <- NULL
}

if (!is.null(weight.censor) & diagnostic!=2)  {
  list.weight.censor   <- paste(weight.censor,times.exposure,sep="_")
} else {list.weight.censor <- NULL
}

if (!is.null(censor) & diagnostic==2)  {
  list.censor   <- paste(censor,times.covariate,sep="_")
} else {list.censor <- NULL
}

if (!is.null(weight.censor) & diagnostic==2)  {
  list.weight.censor   <- paste(weight.censor,times.covariate,sep="_")
} else {list.weight.censor <- NULL
}

if (!is.null(strata))  {
  list.strata   <- paste(strata,times.exposure,sep="_")
} else {list.strata <- NULL
}

if (censoring=="no") {
  censor.unique <- NULL
} else if (censoring=="yes") {
  censor.unique <- censor
}

if (any(!list.exposure %in% names(input))) {
  stop("ERROR: The exposure root name is misspelled, or some exposure measurements are missing from the input dataframe, or incorrect measurement times have been specified")
}

if (!is.null(history) & any(!list.history %in% names(input))) {
  stop("ERROR: Either the history root name is misspelled or some history measurements are missing from the input dataframe.")
}

if (!is.null(censor) & any(!list.censor %in% names(input))) {
  stop("ERROR: Either the censor root name is misspelled or some censor measurements are missing from the input dataframe.")
}

if (!is.null(weight.censor) & any(!list.weight.censor %in% names(input))) {
  stop("ERROR: Either the weight.censor root name is misspelled or some weight.censor measurements are missing from the input dataframe.")
}

if (!is.null(weight.exposure) & any(!list.weight.exposure %in% names(input))) {
  stop("ERROR: Either the weight.exposure root name is misspelled or some weight.exposure measurements are missing from the input dataframe.")
}


  covariate.unique <- c(static.covariate,temporal.covariate)

  #issue a warning if delimiter is contained in root names
  if (any(grepl("_",covariate.unique))) {
  stop("Please ensure that covariate root names do not contain an underscore i.e. '_' ")
  } 
  
  if (is.null(static.covariate)) {
    list.static.covariate <- NULL
	} else if (!is.null(static.covariate)) {
    list.static.covariate <- sort(as.vector(sapply(static.covariate,paste,min(times.covariate,times.exposure),sep="_")))
	}

  if (is.null(temporal.covariate)) {
    list.temporal.covariate <- NULL
	} else if (!is.null(temporal.covariate)) {
    list.temporal.covariate <- sort(as.vector(sapply(temporal.covariate,paste,times.covariate,sep="_")))
	}

  list.covariate <- c(list.static.covariate,list.temporal.covariate)
  list.all.covariate <- names(input)
  list.covariate <- intersect(list.covariate,list.all.covariate)  #remove missing covariate measurements from list
  
  #issue a warning if the exposure and/or covariates contain missing data

  list.exposure.check <- paste(list.exposure, collapse="|")
  list.covariate.check <- paste(list.covariate, collapse="|")

  expCheck <- input[,grep(list.exposure.check, names(input), value=TRUE)]
  covCheck <- input[,grep(list.covariate.check, names(input), value=TRUE)]

  if ( any(is.na(expCheck)) | any(is.na(covCheck))) {
    warning("The exposure and/or some covariates contain missing data. Subsequent calculations are not guaranteed to be unbiased in the presence of partially missing data.")
  }


  #issue an error and abort the program if any exposure or covariate is not numeric

  expCovList <- c(list.exposure, list.covariate)
  covFormat_check <- input[expCovList]

  if ( any(!sapply(covFormat_check,class) %in% "integer") & any(!sapply(covFormat_check,class) %in% "numeric")) {
    stop("ERROR: At least one exposure or covariate is not formatted properly. Please ensure that these variables are in numeric format.")
  }

step1 <- input[,c(id,list.exposure,list.covariate,list.history,list.weight.exposure,list.weight.censor,list.strata,list.censor)]

if (censoring=="no" | (censoring=="yes" & diagnostic!=2)) {
  step2 <- step1 %>% gather_(key="wide.name.exp",value="value.exp",c(list.exposure,list.history,list.weight.exposure,list.weight.censor,list.censor,list.strata))
} else if (censoring=="yes" & diagnostic==2) {
  step2 <- step1 %>% gather_(key="wide.name.exp",value="value.exp",c(list.exposure,list.history,list.weight.exposure,list.censor,list.strata))
}

step3 <- step2 %>% separate(wide.name.exp,c("name.exp","time.exposure"),sep="_") %>%
  spread(key="name.exp",value="value.exp")

if (censoring=="no" | (censoring=="yes" & diagnostic!=2)) {
  step4 <- step3 %>% gather_(key="wide.name.cov",value="value.cov",list.covariate) %>%
    separate(wide.name.cov,c("name.cov","time.covariate"),sep="_")
} else if (censoring=="yes" & diagnostic==2) {
  step4 <- step3 %>% gather_(key="wide.name.cov",value="value.cov",c(list.covariate,list.weight.censor)) %>%
    separate(wide.name.cov,c("name.cov","time.covariate"),sep="_") %>%
    spread(key="name.cov",value="value.cov") %>%
    gather_(key="name.cov",value="value.cov",covariate.unique)
}

if (is.null(list.censor)) {
  step5 <- step4
} else {
  step5 <- step4 %>% rename_("censor"=censor.unique) %>%
    filter(censor==0)
}

#issue warning if any of the inputted covariates are not in the final name.cov column
if ( any(!covariate.unique %in% as.vector(unique(step5$name.cov)))) warning("Some covariates listed in temporal.covariate and/or static.covariate are either absent or have missing values at all timepoints. These covariates will not be included in the balance table or plot.")  

step5 %>% mutate(time.exposure=as.numeric(time.exposure),time.covariate=as.numeric(time.covariate)) %>%
  na.omit()
  
}


################
##OMIT HISTORY##
################

#function to remove irrelevant covariate history from balance tables and plots
omit.history <- function (input,
                          omission,
                          covariate.name,
                          distance=NULL,
                          times=NULL) {

  	
						  
  if (class(input$name.cov) %in% "factor") {		
  check <- "is.factor"
  sort.order <- levels(input$name.cov)  
  input$name.cov <- as.character(input$name.cov)					  
  } else if (class(input$name.cov) %in% "character"){
  check <- "is.character"
  }
  
  if (is.null(omission)) {
    stop("ERROR: 'omission' needs to be specified")
  }

  if (omission=="fixed" && is.null(times)) {
    stop("ERROR: 'times' needs to be specified")
  } else if (omission=="fixed" && any(!times %in% input$time.covariate)) {
    stop("ERROR: one or more values in the 'times' vector are not covariate measurement times in the dataframe")
  }

  if (omission=="relative" && is.null(distance)) {
    stop ("ERROR: 'distance' needs to be specified")
  } else if (omission=="relative" && all(!distance %in% (input$time.exposure-input$time.covariate)))   {
    stop ("ERROR: 'distance' does not equal any difference between time.exposure minus time.covariate in the dataframe.")
  }

  if (omission=="relative") {
    output <- mutate(input,name.cov=ifelse(name.cov %in% covariate.name & (time.exposure-time.covariate)>=distance,NA,name.cov))
  } else if (omission=="fixed") {
    output <- mutate(input,name.cov=ifelse(name.cov %in% covariate.name & (time.covariate %in% times),NA,name.cov))
  } else if (omission=="same.time") {
    output <- mutate(input,name.cov=ifelse(name.cov %in% covariate.name & (time.exposure==time.covariate),NA,name.cov))
  }
 
  if (check=="is.factor") { 
  output <- mutate(output,name.cov=factor(output$name.cov,levels=sort.order))
  }
  
  output <- output %>% na.omit()
}

######################
#APPLY.SCOPE FUNCTION#
######################

apply.scope <- function (	input,
							diagnostic,
							approach,
							scope="all",
							average.over=NULL,
							periods=NULL,
							list.distance=NULL,
							recency=NULL,
							sort.order="alphabetical",
							ignore.missing.metric="no",		
							metric="SMD"               		
							) {     	
	
	if (ignore.missing.metric=="yes") {						
	
	  if (metric=="SMD") {       								
	  input <- input %>% filter(!is.na(SMD))                	
	  } else if (metric=="D") {  								
	  input <- input %>% filter(!is.na(D))                  	
	  }
	}  else if (ignore.missing.metric=="no") {
	}
	
	  if (scope=="all") {

		final.table <- input

	  } else if (scope=="average") {

		if (average.over=="values" | average.over=="strata" | average.over=="history" | average.over=="time" | average.over=="distance" | is.null(average.over)) {

		  if (length(unique(input$E))==1) {
			final.table <- input
		  } else if (length(unique(input$E))!=1) {

			if (approach!="stratify") {
			  grouped.table <- input %>% group_by(H,time.exposure,time.covariate,name.cov)
			} else if (approach=="stratify" & diagnostic==2) {
			  grouped.table <- input %>% group_by(S,time.exposure,time.covariate,name.cov)
			} else if (approach=="stratify" & diagnostic==3) {
			  grouped.table <- input %>% group_by(S,H,time.exposure,time.covariate,name.cov)
			}

			final.table <- grouped.table %>% summarise(D=sum(D*Nexp)/sum(Nexp),
													   SMD=sum(SMD*Nexp)/sum(Nexp),
													   Nexp=sum(Nexp),
													   N=sum(N))
		  }
		}

		if (average.over=="strata" | average.over=="history" | average.over=="time" | average.over=="distance" | is.null(average.over)) {

		  if ("S" %in% names(final.table)) {

			if (diagnostic==2) {
			  grouped.table <- final.table %>% group_by(time.exposure,time.covariate,name.cov)
			} else if (diagnostic==3) {
			  grouped.table <- final.table %>% group_by(H,time.exposure,time.covariate,name.cov)
			}

			final.table <- grouped.table %>% summarise(D=sum(D*N)/sum(N),
													   SMD=sum(SMD*N)/sum(N),
													   Nexp=sum(Nexp),
													   N=sum(N))

		  }
		}

		if (average.over=="history" | average.over=="time" | average.over=="distance" | is.null(average.over)) {

		  if ("H" %in% names(final.table)) {

			final.table <- final.table %>% group_by(time.exposure,time.covariate,name.cov) %>%
			  summarise(D=sum(D*N)/sum(N),
						SMD=sum(SMD*N)/sum(N),
						N=sum(N))

		  }
		}

		if (average.over=="time" | average.over=="distance" | is.null(average.over)) {

		  if (diagnostic==1 | diagnostic==3) {

			final.table <- mutate(final.table,distance=time.exposure-time.covariate,time=time.covariate)

		  } else if (diagnostic==2) {

			final.table <- mutate(final.table,distance=time.covariate-time.exposure,time=time.exposure)

		  }

		  final.table <- final.table %>% group_by(distance,name.cov) %>%
			summarise(D=sum(D*N)/sum(N),
					  SMD=sum(SMD*N)/sum(N),
					  N=sum(N))

		  if (!is.null(list.distance)) {
		  final.table <- final.table %>% filter (distance %in% list.distance)
		  }
		}

		if (average.over=="distance") {

		  make.period <- function (x,y) {

			x <- as.matrix(x)
			z <- matrix(NA,nrow=nrow(x),ncol=3)
			for (i in 1:nrow(x)) {
			  for (j in 1:length(y)) {
				if (x[i,] %in% y[[j]]) {
				  z[i,1] <- j
				  z[i,2] <- min(y[[j]])
				  z[i,3] <- max(y[[j]])
				  z <- data.frame(z)
				}
			  }
			}
			colnames(z) <- c("period.id","period.start","period.end")
			return(z)
		  }

		  if (is.null(periods)) {
			periods[[1]] <- unique(final.table$distance)
		  }
		  
		  period.table <- final.table %>% ungroup() %>% select(distance) %>% make.period(periods)
		  
		  final.table <- data.frame(period.table,final.table) %>% group_by(period.id) %>% 
				mutate(period.start=min(distance),period.end=max(distance)) %>%           
					ungroup()															  
		  
		  final.table <- final.table %>% group_by(period.id,period.start,period.end,name.cov) %>%
			summarise(D=sum(D*N)/sum(N),
					  SMD=sum(SMD*N)/sum(N),
					  N=sum(N))
		}

	  } else if (scope=="recent") {

		if (is.null(recency)) {

		  if (diagnostic!=2) {
			k <- 0
		  } else if (diagnostic==2) {
			k <- 1
		  }

		} else {
		  k <- recency
		}

		if (diagnostic!=2 ) {

		  final.table <- input %>% filter((time.exposure-k)==time.covariate)

		} else if (diagnostic==2) {

		  final.table <- input %>% filter(time.exposure==(time.covariate-k))

		}

	  }

	  if ((length(sort.order)==1) | is.null(sort.order)) {
	    if (class(final.table$name.cov) %in% "factor") { 
		  sort.order <- levels(final.table$name.cov) 
		  final.table <- final.table %>% ungroup() %>% mutate(name.cov=factor(name.cov,levels=sort.order)) 
		} else if (sort.order=="alphabetical") {
		  sorted.cov.names <- sort(unique(final.table$name.cov))
		  final.table <- final.table %>% ungroup() %>% mutate(name.cov=factor(name.cov,levels=sorted.cov.names))
		}		
	  } else if (length(sort.order)>1 & any(sort.order %in% unique(final.table$name.cov))) {
	  final.table <- final.table %>% ungroup() %>% mutate(name.cov=factor(name.cov,levels=sort.order))
	  }

	  
	  if ("E" %in% names(final.table) & "S" %in% names(final.table) & "H" %in% names(final.table)) {
		final.table <- final.table %>% arrange(E,S,H,name.cov,time.exposure,time.covariate)
		} else if ("E" %in% names(final.table) & "S" %in% names(final.table) & !"H" %in% names(final.table)) {
		final.table <- final.table %>% arrange(E,S,name.cov,time.exposure,time.covariate)
		} else if ("E" %in% names(final.table) & !"S" %in% names(final.table) & "H" %in% names(final.table)) {
		final.table <- final.table %>% arrange(E,H,name.cov,time.exposure,time.covariate)
		} else if (!"E" %in% names(final.table) & "S" %in% names(final.table) & !"H" %in% names(final.table)) {
		final.table <- final.table %>% arrange(S,name.cov,time.exposure,time.covariate)	 
		} else if (!"E" %in% names(final.table) & !"S" %in% names(final.table) & "H" %in% names(final.table)) {
		final.table <- final.table %>% arrange(H,name.cov,time.exposure,time.covariate)	 
		} else if (!"E" %in% names(final.table) & !"S" %in% names(final.table) 
		                                        & !"H" %in% names(final.table) 
												& !"distance" %in% names(final.table)  
												& !"period.id" %in% names(final.table)) {  
		final.table <- final.table %>% arrange(name.cov,time.exposure,time.covariate)															
		} else if ("distance" %in% names(final.table)) {
		final.table <- final.table %>% arrange(distance,name.cov)
		} else if ("period.id" %in% names(final.table)) {
		final.table <- final.table %>% arrange(period.id,period.start,period.end,name.cov)
		}
}

####################
##BALANCE FUNCTION##
####################

#function to create a balance table for a specified diagnostic
balance <- function (input,
                     diagnostic,
                     approach="none",
                     censoring,
                     scope,
                     times.exposure,
                     times.covariate,
                     exposure,
                     history=NULL,
                     weight.exposure=NULL,
                     weight.censor=NULL,
                     strata=NULL,
                     recency=NULL,
                     average.over=NULL,
                     periods=NULL,
					 list.distance=NULL,
					 sort.order="alphabetical",
					 loop="no",
					 ignore.missing.metric="no",   
					 metric="SMD") {               

  if(is.null(input)) {
    stop ("ERROR: 'input' is missing. Please specify the dataframe created by the lengthen() function")
  }
  if(is.null(diagnostic) | !diagnostic %in% c(1,2,3)) {
    stop ("ERROR: 'diagnostic' is missing or misspecified. Please specify as 1, 2 or 3")
  }
  if(is.null(approach) | !approach %in% c("none","weight","stratify")) {
    stop ("ERROR: 'approach' is missing or misspecified. Please specify as none, weight, or stratify")
 }

  if (diagnostic!=1 & approach=="none") stop ("'diagnostic' has been specified as 2 or 3, and to implement the choice properly, approach needs to be specified as weight or stratify.")

  if(is.null(censoring) | !censoring %in% c("no","yes")) {
    stop ("ERROR: 'censoring' is missing. Please specify it as yes or no")
  }
  if (is.null(scope) | !scope %in% c("all","recent","average")) {
    stop ("ERROR: 'scope' is missing. Please specify either all, recent, or average")
  }
  if(is.null(exposure) | !exposure %in% names(input)) {
    stop ("ERROR: Either 'exposure' has not been specified OR the exposure root name is not present in the input dataframe")
  }
  if(!is.null(history) && !history %in% names(input)){
    stop("ERROR: the specified root name for 'history' is not present in the input dataframe")
  }
  if(!is.null(strata) && !strata %in% names(input)){
    stop("ERROR: the specified root name for 'strata' is not present in the input dataframe")
  }
  if(!is.null(weight.exposure) && !weight.exposure %in% names(input)){
    stop("ERROR: the specified root name for 'weight.exposure' is not present in the input dataframe")
  }
  if(!is.null(weight.censor) && !weight.censor %in% names(input)){
    stop("ERROR: the specified root name for 'weight.censor' is not present in the input dataframe")
  }
  if(is.null(times.exposure)) {
    stop ("ERROR: 'times.exposure' is missing. Please specify an integer for times.exposure or a numeric vector of times")
  }
  if(is.null(times.covariate)) {
    stop ("ERROR: 'times.covariate' is missing. Please specify an integer for times.covariate or a numeric vector of times")
  }
  if (scope=="recent" & is.null(recency)) {
    stop ("ERROR: 'recency' is missing. Please specify an integer between 0 and t for diagnostics 1 and 3, or 0 and t-1 for diagnostic 2")
  }
  if (scope=="average" & (is.null(average.over))) {
    stop ("ERROR: 'average.over' is missing. Please specify one of the following: values, strata, history, time, distance.")
  }
  
  if (!is.null(periods)) {
    
	vector.periods <- c(unlist(periods))
	
	if (!is.list(periods)) { 
	stop ("ERROR: When specifying 'periods', it must be specified as a list e.g. list(0,2:3,5:9)") 
	} else if ((!class(vector.periods) %in% "numeric") & (!class(vector.periods) %in% "integer")) { 
    stop ("ERROR: When specifying 'periods', the list must only contain numeric or integer values e.g. list(0,2:3,5:9)") 
	} else if (!all(!vector.periods %in% vector.periods[duplicated(vector.periods)])) {
	stop ("ERROR: When specifying 'periods', the values should be unique e.g. list(0,2:3,5:9") 
	} else if ((!any(vector.periods)==sort(vector.periods))) { 
	stop ("ERROR: When specifying 'periods', the values must be in sorted order e.g. list(0,2:3,5:9)") 
	}
  }
	
  if (diagnostic==1) {
    if (is.null(history)) stop ("ERROR: 'history' is missing. Please specify a root name for exposure history")
  }  else if ((diagnostic==2 | diagnostic==3) & ((is.null(history) | is.null(weight.exposure)) & is.null(strata))) {
    stop ("ERROR: Please specify the root names for exposure history and exposure weight variables, or the root name for strata variables")
  } else if ((diagnostic==3 & approach=="stratify") & (is.null(history) | is.null(strata))) {
  stop ("ERROR: Please specify the root names for exposure history and/or strata variables")
  }

  if ((length(sort.order)==1) | is.null(sort.order)) {
    if (sort.order!="alphabetical") {
    stop ("ERROR: either specify sort.order as 'alphabetical' or as a character vector of covariates")
    }
  } else if (length(sort.order)>1 & !all(unique(input$name.cov) %in% unique(sort.order))) {
  stop ("ERROR: when specifying the character vector for sort.order, include all covariate names in the input dataframe, and also ensure that their spelling match those specified in sort.order. Provided these criteria are met, the software will still run even if sort.order includes extraneous covariate names that are NOT present in the input dataframe")
  } 
  
  if (is.null(history)) {
    input$history.none <- "H"
    history <- "history.none"
  }

  if (is.null(strata)) {
    input$strata.none <- 1
    strata <- "strata.none"
  }

  if (is.null(weight.exposure)) {
    input$weight.exposure.none <- 1
    weight.exposure <- "weight.exposure.none"
  }

  if (is.null(weight.censor)) {
    input$weight.censor.none <- 1
    weight.censor <- "weight.censor.none"
  }

  t.exp.data <- unique(input$time.exposure)
  t.cov.data <- unique(input$time.covariate)
  t.exp.spec <- unique(times.exposure)
  t.cov.spec <- unique(times.covariate)

  if ((all(t.exp.data %in% t.exp.spec)) & (all(t.cov.data %in% t.cov.spec))) {
  } else {
    input <- input %>% filter((time.exposure %in% t.exp.spec) & (time.covariate %in% t.cov.spec))
  }

  input <- rename_(input,"E"=exposure,"H"=history,"S"=strata,"W_a"=weight.exposure,"W_s"=weight.censor)

  if (approach=="weight" | approach=="none") {

    if (censoring=="yes") {
      input <- mutate(input,W=as.numeric(W_a)*as.numeric(W_s))
    } else if (censoring=="no") {
      input <- mutate(input,W=as.numeric(W_a))
    }

    if (diagnostic==1) {

      input <- mutate(input,W=1)

      temp.table <-
        data.frame(input %>% select(E,H,W,time.exposure,time.covariate,name.cov,value.cov) %>%
                     group_by(E,H,name.cov,time.exposure,time.covariate) %>%
                     summarise(mean.cov_b=weighted.mean(x=value.cov,w=W,na.rm=TRUE),
                               sd.cov_b=sd(x=value.cov,na.rm=TRUE),
                               n.cov_b=sum(W)))

        check.table <- temp.table %>%
          group_by(H,time.exposure) %>%
            summarise(nexpval=n_distinct(E)) %>%
              group_by()

        if (all(check.table$nexpval==1)) stop("ERROR: None of the exposure times have exposure variation within levels of exposure history. The program has terminated because the resulting balance table is empty")
        if (any(check.table$nexpval==1)) warning("Some exposure times have no exposure variation within levels of exposure history. Estimates for these times will not appear in the results")

        temp.table <- temp.table %>% group_by(H,name.cov,time.exposure,time.covariate) %>%
                       mutate(mean.cov_a=first(mean.cov_b),
                            sd.cov_a=first(sd.cov_b),
                            n.cov_a=first(n.cov_b)) %>%
                       filter(E!=first(E))

      full.table <- temp.table %>% mutate(D=mean.cov_b-mean.cov_a,
                                          SMD=ifelse(D==0,0,
                                                ifelse(sd.cov_b==0 | sd.cov_a==0,NA_real_,
												(mean.cov_b-mean.cov_a)/sqrt((sd.cov_a^2*(n.cov_a-1)+sd.cov_b^2*(n.cov_b-1))/(n.cov_a+n.cov_b-2)))),
                                          N=n.cov_a+n.cov_b,
                                          Nexp=n.cov_b)

      if ( any(is.na(unique(full.table$SMD))) ) warning("SMD values have been set to missing where there is no covariate variation within some level of time-exposure, time-covariate, exposure history, and exposure value; in this case averages for SMD estimates will also appear as missing")
      if ( all(full.table$D==0) & all(full.table$SMD==0) ) warning("There may be no covariate variation within any level of time-exposure, time-covariate, exposure history and/or strata, and exposure value; please ensure that the temporal covariates are specified correctly.")


      sub.table  <- full.table %>%
                      select(E,H,name.cov,time.exposure,time.covariate,D,SMD,N,Nexp) %>%
                        filter (!is.na(D)) %>%
                          filter(time.exposure>=time.covariate) %>%
                            arrange(name.cov,time.exposure,time.covariate,H)

    } else if (diagnostic==2 | diagnostic==3) {

      temp.table <-
        data.frame(input %>% select(E,H,W,time.exposure,time.covariate,name.cov,value.cov) %>%
                     group_by(E,H,name.cov,time.exposure,time.covariate) %>%
                     summarise(mean.cov_b=weighted.mean(x=value.cov,w=W,na.rm=TRUE),
                               sd.cov_b=sd(x=value.cov,na.rm=TRUE),
                               n.cov_b=sum(W)))

        check.table <- temp.table %>%
          group_by(H,time.exposure) %>%
            summarise(nexpval=n_distinct(E)) %>%
              group_by()

        if (all(check.table$nexpval==1)) stop("ERROR: None of the exposure times have exposure variation within levels of exposure history. The program has terminated because the resulting balance table is empty")
        if (any(check.table$nexpval==1)) warning("Some exposure times have no exposure variation within levels of exposure history. Estimates for these times will not appear in the results")

        temp.table <- temp.table %>% group_by(H,name.cov,time.exposure,time.covariate) %>%
                     mutate(mean.cov_a=first(mean.cov_b),
                            sd.cov_a=first(sd.cov_b),
                            n.cov_a=first(n.cov_b)) %>%
                     filter(E!=first(E))

      full.table <- temp.table %>% mutate(D=mean.cov_b-mean.cov_a,
                                          SMD=ifelse(D==0,0,
                                                     ifelse(sd.cov_b==0 | sd.cov_a==0,NA_real_,
                                                            (mean.cov_b-mean.cov_a)/sqrt((sd.cov_a^2*(n.cov_a-1)+sd.cov_b^2*(n.cov_b-1))/(n.cov_a+n.cov_b-2)))),
                                          N=n.cov_a+n.cov_b,
                                          Nexp=n.cov_b)

      if ( any(is.na(unique(full.table$SMD))) ) warning("SMD values have been set to missing where there is no covariate variation within  some level of time-exposure, time-covariate, exposure history, and exposure value; in this case averages for SMD estimates will also appear as missing")
      if ( all(full.table$D==0) & all(full.table$SMD==0) ) warning("There may be no covariate variation within any level of time-exposure, time-covariate, exposure history and/or strata, and exposure value; please ensure that the temporal covariates are specified correctly.")


      full.table <- full.table %>% select(E,H,name.cov,time.exposure,time.covariate,D,SMD,N,Nexp)

      if (diagnostic==2) {

        sub.table <- full.table %>% filter (!is.na(D)) %>% filter(time.exposure<time.covariate)

      } else if (diagnostic==3) {

        sub.table <- full.table %>% filter (!is.na(D)) %>% filter(time.exposure>=time.covariate)
      }
    }

  } else if (approach=="stratify") {

    if (censoring=="yes") {
      input <- mutate(input,W=as.numeric(W_s))
    } else if (censoring=="no") {
      input <- mutate(input,W=1)
    }

    if (diagnostic==1) {

      input <- mutate(input,W=1)

      temp.table <-
        data.frame(input %>% select(E,H,W,time.exposure,time.covariate,name.cov,value.cov) %>%
                     group_by(E,H,name.cov,time.exposure,time.covariate) %>%
                     summarise(mean.cov_b=weighted.mean(x=value.cov,w=W,na.rm=TRUE),
                               sd.cov_b=sd(x=value.cov,na.rm=TRUE),
                               n.cov_b=sum(W)))

        check.table <- temp.table %>%
          group_by(H,time.exposure) %>%
            summarise(nexpval=n_distinct(E)) %>%
              group_by()

        if (all(check.table$nexpval==1)) stop("ERROR: None of the exposure times have exposure variation within levels of exposure history. The program has terminated because the resulting balance table is empty")
        if (any(check.table$nexpval==1)) warning("Some exposure times have no exposure variation within levels of exposure history. Estimates for these times will not appear in the results")

        temp.table <- temp.table %>% group_by(H,name.cov,time.exposure,time.covariate) %>%
                     mutate(mean.cov_a=first(mean.cov_b),
                            sd.cov_a=first(sd.cov_b),
                            n.cov_a=first(n.cov_b)) %>%
                     filter(E!=first(E))

      full.table <- temp.table %>% mutate(D=mean.cov_b-mean.cov_a,
                                          SMD=ifelse(D==0,0,
                                                     ifelse(sd.cov_b==0 | sd.cov_a==0,NA_real_,
                                                            (mean.cov_b-mean.cov_a)/sqrt((sd.cov_a^2*(n.cov_a-1)+sd.cov_b^2*(n.cov_b-1))/(n.cov_a+n.cov_b-2)))),
                                          N=n.cov_a+n.cov_b,
                                          Nexp=n.cov_b)

      if ( any(is.na(unique(full.table$SMD))) ) warning("SMD values have been set to missing where there is no covariate variation within some level of time-exposure, time-covariate, exposure history, and exposure value; in this case averages for SMD estimates will also appear as missing")
      if ( all(full.table$D==0) & all(full.table$SMD==0) ) warning("There may be no covariate variation within any level of time-exposure, time-covariate, exposure history and/or strata, and exposure value; please ensure that the temporal covariates are specified correctly.")


      sub.table <- full.table %>%
                    select(E,H,name.cov,time.exposure,time.covariate,D,SMD,N,Nexp) %>%
                      filter (!is.na(D)) %>%
                        filter(time.exposure>=time.covariate) %>%
                          arrange(name.cov,time.exposure,time.covariate,H)

    } else if (diagnostic==2) {

      values.exposure <- sort(unique(input$E))
      temp.table <-
        data.frame(input %>% select(E,S,W,time.exposure,time.covariate,name.cov,value.cov) %>%
                     group_by(E,S,name.cov,time.exposure,time.covariate) %>%
                     summarise(mean.cov_b=weighted.mean(x=value.cov,w=W,na.rm=TRUE),
                               sd.cov_b=sd(x=value.cov,na.rm=TRUE),
                               n.cov_b=sum(W)))

       check.table <- temp.table %>%
         group_by(S,time.exposure) %>%
           summarise(nexpval=n_distinct(E)) %>%
             group_by()

       if (all(check.table$nexpval==1)) stop("ERROR: None of the exposure times have exposure variation within levels of exposure history. The program has terminated because the resulting balance table is empty")
       if (any(check.table$nexpval==1)) warning("Some exposure times have no exposure variation within levels of strata. Estimates for these times will not appear in the results")

        temp.table <- temp.table %>% group_by(S,name.cov,time.exposure,time.covariate) %>%
                     mutate(mean.cov_a=first(mean.cov_b),
                            sd.cov_a=first(sd.cov_b),
                            n.cov_a=first(n.cov_b)) %>%
                     filter(E!=first(E))

      full.table <- temp.table %>% mutate(D=mean.cov_b-mean.cov_a,
                                          SMD=ifelse(D==0,0,
                                                     ifelse(sd.cov_b==0 | sd.cov_a==0,NA_real_,
                                                            (mean.cov_b-mean.cov_a)/sqrt((sd.cov_a^2*(n.cov_a-1)+sd.cov_b^2*(n.cov_b-1))/(n.cov_a+n.cov_b-2)))),
                                          N=n.cov_a+n.cov_b,
                                          Nexp=n.cov_b)

      if ( any(is.na(unique(full.table$SMD))) ) warning("SMD values have been set to missing where there is no covariate variation within some level of time-exposure, time-covariate, strata, and exposure value; in this case averages for SMD estimates will also appear as missing")
      if ( all(full.table$D==0) & all(full.table$SMD==0) ) warning("There may be no covariate variation within any level of time-exposure, time-covariate, exposure history and/or strata, and exposure value; please ensure that the temporal covariates are specified correctly.")


      sub.table <-  full.table %>%
                     select(E,S,name.cov,time.exposure,time.covariate,D,SMD,N,Nexp) %>%
                       filter (!is.na(D)) %>%
                         filter(time.exposure<time.covariate) %>%
                           arrange(S,name.cov,time.exposure,time.covariate)

    } else if (diagnostic==3) {

      temp.table <-
        data.frame(input %>% select(E,S,H,W,time.exposure,time.covariate,name.cov,value.cov) %>%
                     group_by(E,S,H,time.exposure,time.covariate,name.cov) %>%
                     summarise(mean.cov_b=weighted.mean(x=value.cov,w=W,na.rm=TRUE),
                               sd.cov_b=sd(x=value.cov,na.rm=TRUE),
                               n.cov_b=sum(W)))

        check.table <- temp.table %>%
          group_by(S,H,time.exposure) %>%
            summarise(nexpval=n_distinct(E)) %>%
              group_by()

        if (all(check.table$nexpval==1)) stop("ERROR: None of the exposure times have exposure variation within levels of exposure history. The program has terminated because the resulting balance table is empty")
        if (any(check.table$nexpval==1)) warning("Some exposure times have no exposure variation within levels of strata and exposure history. Estimates for these times will not appear in the results")

        temp.table <- temp.table %>% group_by(S,H,name.cov,time.exposure,time.covariate) %>%
                     mutate(mean.cov_a=first(mean.cov_b),
                            sd.cov_a=first(sd.cov_b),
                            n.cov_a=first(n.cov_b)) %>%
                     filter(E!=first(E))

       full.table <- temp.table %>% mutate(D=mean.cov_b-mean.cov_a,
                                          SMD=ifelse(D==0,0,
                                                     ifelse(sd.cov_b==0 | sd.cov_a==0,NA_real_,
                                                            (mean.cov_b-mean.cov_a)/sqrt((sd.cov_a^2*(n.cov_a-1)+sd.cov_b^2*(n.cov_b-1))/(n.cov_a+n.cov_b-2)))),
                                          N=n.cov_a+n.cov_b,
                                          Nexp=n.cov_b)

        if ( any(is.na(unique(full.table$SMD))) ) warning("SMD values have been set to missing where there is no covariate variation within some level of time-exposure, time-covariate, strata, exposure history, and exposure value; in this case averages for SMD estimates will also appear as missing")
        if ( all(full.table$D==0) & all(full.table$SMD==0) ) warning("There may be no covariate variation within any level of time-exposure, time-covariate, exposure history and/or strata, and exposure value; please ensure that the temporal covariates are specified correctly.")


        sub.table <- full.table %>%
                       select(E,S,H,name.cov,time.exposure,time.covariate,D,SMD,N,Nexp) %>%
                         filter (!is.na(D)) %>%
                           filter(time.exposure>=time.covariate) %>%
                             arrange(S,name.cov,time.exposure,time.covariate,H)

    }
  }
  
  if (loop=="no") {  
  output <- apply.scope(input=sub.table, 
						diagnostic=diagnostic, 
						approach=approach, 
						scope=scope, 
						average.over=average.over, 
						periods=periods, 
						list.distance=list.distance, 
						recency=recency, 
						sort.order=sort.order, 
						ignore.missing.metric=ignore.missing.metric, 
						metric=metric) 
  } else if (loop=="yes") {  
  output <- sub.table
  }
  return(output)
}



##########
#DIAGNOSE#
##########

diagnose <- function (
				input,
				diagnostic,
				approach="none",
				scope,
				censoring,
				id,
				times.exposure,
				times.covariate,
				exposure,
				temporal.covariate,
				static.covariate=NULL,
				history=NULL,
				weight.exposure=NULL,
				censor=NULL,
				weight.censor=NULL,
				strata=NULL,
				recency=NULL,
				average.over=NULL,
				periods=NULL,
				list.distance=NULL,
				sort.order="alphabetical",
				loop="no",
				loop.type=NULL,
				ignore.missing.metric="no",
				metric="SMD"
) {

	loop.fxn <- function (arg.temporal.covariate,arg.static.covariate,arg.times.exposure,arg.times.covariate) {

		output.lengthen <- lengthen (
			input=input,
			diagnostic=diagnostic,
			censoring=censoring,
			id=id,
			times.exposure=arg.times.exposure,
			times.covariate=arg.times.covariate,
			exposure=exposure,
			temporal.covariate=arg.temporal.covariate,
			static.covariate=arg.static.covariate,
			history=history,
			weight.exposure=weight.exposure,
			censor=censor,
			weight.censor=weight.censor,
			strata=strata
			)

		output.balance <- balance (
			input=output.lengthen,
			diagnostic=diagnostic,
			approach=approach,
			censoring=censoring,
			scope=scope,
			times.exposure=arg.times.exposure,
			times.covariate=arg.times.covariate,
			exposure=exposure,
			history=history,
			weight.exposure=weight.exposure,
			weight.censor=weight.censor,
			strata=strata,
			recency=recency,
			average.over=average.over,
			periods=periods,
			list.distance=list.distance,
			sort.order=sort.order,
			loop=loop,
			ignore.missing.metric=ignore.missing.metric,	
			metric=metric								   
			)

	}


	if (loop=="no") {

				output <- loop.fxn (
							arg.temporal.covariate=temporal.covariate,
							arg.static.covariate=static.covariate,
							arg.times.exposure=times.exposure,
							arg.times.covariate=times.covariate
							)
							

	} else if (loop=="yes") {
	
	    if (is.null(loop.type)) {
		stop("Please specify the type of loop: 1 to loop over covariates; 2 to loop over covariates and exposure times; 3 to loop over covariates and exposure/covariate measurement time pairs")
		}
	
		if (loop.type==1) { #LOOP OVER COVARIATES ONLY
		
			results.temporal <- list()
			
			for (I in seq_along(temporal.covariate)) {
				results.temporal[[I]] <- loop.fxn (
							arg.temporal.covariate=temporal.covariate[I],
							arg.static.covariate=NULL,
							arg.times.exposure=times.exposure,
							arg.times.covariate=times.covariate
							)			

			}	

			results.static <- list()
			
			for (L in seq_along(static.covariate)) {
				results.static[[L]] <- loop.fxn (
							arg.temporal.covariate=NULL,
							arg.static.covariate=static.covariate[L],
							arg.times.exposure=times.exposure,
							arg.times.covariate=times.covariate
							)			
			}		
		
			results.all <- c(results.temporal,results.static)
			output <- results.all %>% 
				rbindlist() %>% 
					as.data.frame() %>%
						apply.scope (
							diagnostic=diagnostic,
							approach=approach,
							scope=scope,
							average.over=average.over,
							periods=periods,
							list.distance=list.distance,
							recency=recency,
							sort.order=sort.order,
							ignore.missing.metric=ignore.missing.metric,
							metric=metric
							)	
						
		} else if (loop.type==2) { #LOOP OVER COVARIATES AND EXPOSURE TIMES (DIAGNOSTIC 1|3) OR COVARIATE TIMES (DIAGNOSTIC 2)
		
			if (diagnostic==1 | diagnostic==3) {
			loop.order <- "exogeneity"
			loop.sequence.j <- times.exposure
			} else if (diagnostic==2) {
			loop.order <- "feedback"
			loop.sequence.j <- times.covariate
			}		

			results.temporal <- list()
			results.temporal.iter <- list()
			
			for (I in seq_along(temporal.covariate)) {
				
				for (J in seq_along(loop.sequence.j)) {			
			
					if (loop.order=="exogeneity") {
					iter.times.exposure  <- loop.sequence.j[J]
					iter.times.covariate <- times.covariate[which(times.covariate <= iter.times.exposure)]
					} else if (loop.order=="feedback") {
					iter.times.covariate <- loop.sequence.j[J]
					iter.times.exposure  <- times.exposure[which(times.exposure < iter.times.covariate)]				
					}
			
					results.temporal.iter[[J]] <- loop.fxn (
							arg.temporal.covariate=temporal.covariate[I],
							arg.static.covariate=NULL,
							arg.times.exposure=iter.times.exposure,
							arg.times.covariate=iter.times.covariate
							)
								
				}
				results.temporal <- c(results.temporal,results.temporal.iter)
			}	
			
			if (loop.order=="exogeneity") {

				results.static <- list()
				results.static.iter <- list()
				
				min.times.covariate <- min(times.covariate)
				
				for (L in seq_along(static.covariate)) {
				
					for (M in seq_along(times.exposure)) {
				
					results.static.iter[[M]] <- loop.fxn (
								arg.temporal.covariate=NULL,
								arg.static.covariate=static.covariate[L],
								arg.times.exposure=times.exposure[M],
								arg.times.covariate=min.times.covariate
								)			
					}
					results.static <- c(results.static,results.static.iter)			
				}		
			} else if (loop.order=="feedback") {
			results.static <- NULL
			}
			
			results.all <- c(results.temporal,results.static)
			output <- results.all %>% 
				rbindlist() %>% 
					as.data.frame() %>%
						apply.scope (
							diagnostic=diagnostic,
							approach=approach,
							scope=scope,
							average.over=average.over,
							periods=periods,
							list.distance=list.distance,
							recency=recency,
							sort.order=sort.order,
							ignore.missing.metric=ignore.missing.metric,
							metric=metric
							)		
		
		} else if (loop.type==3) { #LOOP OVER COVARIATES, EXPOSURE TIMES AND COVARIATE TIMES

			if (diagnostic==1 | diagnostic==3) {
			loop.order <- "exogeneity"
			loop.sequence.j <- times.exposure
			} else if (diagnostic==2) {
			loop.order <- "feedback"
			loop.sequence.j <- times.covariate
			}
		
			results.temporal <- list()
			results.temporal.iter <- list()
			
			for (I in seq_along(temporal.covariate)) {
				
				for (J in seq_along(times.exposure)) {	

					for (K in seq_along(times.covariate)) {
			
						if (loop.type=="exogeneity" & times.exposure[J] < times.covariate[K]) {
						next
						} else if (loop.type=="feedback" & times.exposure[J] >= times.covariate[K]) {
						next
						} else {
			
						results.temporal.iter[[K]] <- loop.fxn (
								arg.temporal.covariate=temporal.covariate[I],
								arg.static.covariate=NULL,
								arg.times.exposure=times.exposure[J],
								arg.times.covariate=times.covariate[K]
								)
						}			
					}
					results.temporal <- c(results.temporal,results.temporal.iter)
				}
			}	
	

			if (loop.order=="exogeneity") {
			
				min.times.covariate <- min(times.covariate)
			
				results.static <- list()
				results.static.iter <- list()
				
				for (L in seq_along(static.covariate)) {
				
					for (M in seq_along(times.exposure)) {
				
					results.static.iter[[M]] <- loop.fxn (
								arg.temporal.covariate=NULL,
								arg.static.covariate=static.covariate[L],
								arg.times.exposure=times.exposure[M],
								arg.times.covariate=min.times.covariate
								)			
					}
					results.static <- c(results.static,results.static.iter)			
				}

			} else if (loop.order=="feedback") {
			results.static <- NULL
			}
				
			results.all <- c(results.temporal,results.static)
			output <- results.all %>% 
				rbindlist() %>% 
					as.data.frame() %>%
						apply.scope (
							diagnostic=diagnostic,
							approach=approach,
							scope=scope,
							average.over=average.over,
							periods=periods,
							list.distance=list.distance,
							recency=recency,
							sort.order=sort.order,
							ignore.missing.metric=ignore.missing.metric,
							metric=metric
							)
						
		}

 
	}
	return(output)	
}






#####################
##MAKEPLOT FUNCTION##
#####################

#function to create balance plot for a specified diagnostic
makeplot <- function (input,
                      diagnostic,
                      approach,
                      metric="SMD",
                      scope,
                      average.over=NULL,
                      stratum=NULL,
                      label.exposure="A",
                      label.covariate="C",
                      lbound=-1,
                      ubound=1,
                      ratio=2,
                      text.axis.title=8,
                      text.axis.y=6.5,
                      text.axis.x=6.5,
                      text.strip.y=10,
                      text.strip.x=10,
                      point.size=.75,
                      zeroline.size=.1,
                      refline.size=.1,
                      refline.limit.a=-.25,
                      refline.limit.b=0.25,
                      panel.margin.size=.75,
                      axis.title=NULL,
					  label.width=15) {

  if(is.null(input)) {
    stop ("ERROR: 'input' is missing. Please specify the dataframe created by the balance() function")
  }
  if(is.null(diagnostic) | !diagnostic %in% c(1,2,3)) {
    stop ("ERROR: 'diagnostic' is missing or misspecified. Please specify as 1, 2 or 3")
  }
  if(is.null(approach) | !approach %in% c("none","weight","stratify")) {
    stop ("ERROR: 'approach' is missing or misspecified. Please specify as none, weight, or stratify")
  }
  if(is.null(metric) | !metric %in% c("D","SMD")) {
    stop ("ERROR: 'metric' is missing or misspecified. Please specify as D (for the mean difference) or SMD (for the standardized mean difference)")
  }
  if (is.null(scope) | !scope %in% c("all","recent","average")) {
    stop ("ERROR: 'scope' is missing. Please specify either all, recent, or average")
  }
  if (scope=="average" & (is.null(average.over))) {
    stop ("ERROR: 'average.over' is missing. Please specify one of the following: values, strata, history, time, distance.")
  }
  if (approach=="stratify" & is.null(stratum) & scope!="average") {
    stop ("ERROR: 'stratum' is missing. Please specify an integer indicating the stratum to select for plotting. Note that 'stratum' is not required when strata are averaged over e.g. scope equals average and average.over equals anything higher than values")
  } else if (approach=="stratify" && is.null(stratum) && scope=="average" && (!is.null(average.over) && average.over=="values")) {
    stop ("ERROR: 'stratum' is missing. Please specify a value indicating the stratum to select for plotting. Note that 'stratum' is not required when strata are averaged over e.g. scope equals average and average.over equals anything higher than values")
  } else {
  }

  if (metric=="D") {
    input <- mutate(input,plot.metric=D)
  } else if (metric=="SMD") {
    input <- mutate(input,plot.metric=SMD)
  }

  nonmiss.metric <- sum(!is.na(input$plot.metric))

  if(nonmiss.metric==0) {
    stop ("ERROR: There are no non-missing values for the specified metric. The program has terminated because the resulting plot is empty")
  }


  nonzero.metric <- sum(input$plot.metric!=0)

  if(nonzero.metric==0) {
    warning ("There are no non-zero values for the specified metric. The plot will still be produced, but it is advisable to revisit the balance function and ensure that the temporal covariates were specified correctly")
  }



  if (is.null(stratum) & approach!="stratify") {
    stratum <- rep(1,nrow(input))
  }

  themes <- theme(aspect.ratio=ratio,
                  axis.title = element_text(size=text.axis.title,face="bold"),
                  axis.text.y = element_text(size=text.axis.y,colour="black",vjust=0.33),
                  axis.text.x = element_text(size=text.axis.x,colour="black"),
                  axis.ticks=element_blank(),
                  strip.text.y = element_text(size = text.strip.y),
                  strip.text.x = element_text(size = text.strip.x),
                  strip.background=element_blank(),
                  panel.grid.major.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.y=element_blank(),
                  panel.grid.minor.y=element_blank(),
                  panel.background = element_blank(),
                  panel.margin = unit(panel.margin.size, "lines"))

  if (is.null(average.over)) {
    average.over=""
  }

  if (approach=="weight" | approach=="none" | (approach=="stratify" & scope=="average" & average.over!="values")) {
    sub.input <- input
  } else if (approach=="stratify" & (scope=="all" | scope=="recent" | (scope=="average" & average.over=="values"))) {
    sub.input <- input %>% filter(S==stratum)
  }

  if (scope=="all" | scope=="recent" | (scope=="average" & average.over!="time" & average.over!="distance")) {

    labelled.input <- sub.input %>% mutate(exposure=paste(label.exposure,"(",time.exposure,")",sep=""),
                                           covariate=paste(label.covariate,"(",time.covariate,")",sep=""),
                                           comparison=paste(exposure," vs ",covariate,sep="")
                                           )
	
	values.exposure  <- labelled.input %>% select(time.exposure,exposure) %>% unique()
	values.covariate <- labelled.input %>% select(time.covariate,covariate) %>% unique()
	values.comparison <- labelled.input %>% select(comparison,time.exposure,time.covariate) %>% unique()

	AscendOrderExposure   <- arrange(values.exposure,time.exposure)
	AscendOrderCovariate  <- arrange(values.covariate,time.covariate)
	DescendOrderExposure  <- arrange(values.exposure,desc(time.exposure))
	DescendOrderCovariate <- arrange(values.covariate,desc(time.covariate))
	AscendOrderComparison <- arrange(values.comparison,time.exposure,time.covariate)	
									
    labelled.input <- labelled.input %>% mutate(exposure=factor(exposure,levels=AscendOrderExposure$exposure),
                                                covariate=factor(covariate,levels=AscendOrderCovariate$covariate),
												rev.exposure=factor(exposure,levels=DescendOrderExposure$exposure),
                                                rev.covariate=factor(covariate,levels=DescendOrderCovariate$covariate),
												comparison=factor(comparison,levels=AscendOrderComparison$comparison))

  } else if (average.over=="time" & diagnostic!=2) {
    labelled.input <- sub.input %>% mutate(comparison=paste(label.covariate,"(t-",distance,") vs ",label.exposure,"(t)",sep=""))
	values.distance        <- unique(labelled.input$distance)
	values.comparison      <- labelled.input %>% select(comparison,distance) %>% unique()
	DescendOrderComparison <- arrange(values.comparison,desc(distance),desc(comparison))

  } else if (average.over=="time" & diagnostic==2) {
    labelled.input <- sub.input %>% mutate(comparison=paste(label.exposure,"(t-",distance,") vs ",label.covariate,"(t)",sep=""))
	values.distance        <- unique(labelled.input$distance)
	values.comparison      <- labelled.input %>% select(comparison,distance) %>% unique()
	DescendOrderComparison <- arrange(values.comparison,desc(distance),desc(comparison))

  } else if (average.over=="distance" & diagnostic!=2) {
    labelled.input <- sub.input %>% mutate(comparison=paste(label.covariate,"(t-",period.end,":","t-",period.start,") vs ",label.exposure,"(t)",sep=""),
                                           comparison=ifelse(period.start==period.end,paste(label.covariate,"(t-",period.start,") vs ",label.exposure,"(t)",sep=""),comparison))
	values.period.end      <- unique(labelled.input$period.end)
	values.comparison      <- labelled.input %>% select(period.start,period.end,comparison) %>% unique() 
	DescendOrderComparison <- arrange(values.comparison,desc(period.end),desc(period.start),desc(comparison)) 										   
										   
  } else if (average.over=="distance" & diagnostic==2) {
    labelled.input <- sub.input %>% mutate(comparison=paste(label.exposure,"(t-",period.end,":","t-",period.start,") vs ",label.covariate,"(t)",sep=""),
                                           comparison=ifelse(period.start==period.end,paste(label.exposure,"(t-",period.start,") vs ",label.covariate,"(t)",sep=""),comparison))
	values.period.end      <- unique(labelled.input$period.end)
	values.comparison      <- labelled.input %>% select(period.start,period.end,comparison) %>% unique() 
	DescendOrderComparison <- arrange(values.comparison,desc(period.end),desc(period.start),desc(comparison)) 
  }
  
  if (average.over=="time" | average.over=="distance" | is.null(average.over)) {
    labelled.input <- labelled.input %>% mutate(comparison=factor(comparison,levels=DescendOrderComparison$comparison))
  }


  if ("E" %in% names(labelled.input)) {
    labelled.input <- labelled.input %>% mutate(E=as.factor(E))
  } else {
    labelled.input <- labelled.input %>% mutate(E=as.factor(1))
  }

  if ("H" %in% names(labelled.input)) {
    labelled.input <- labelled.input %>% mutate(H=as.factor(H))
  } else {
    labelled.input <- labelled.input %>% mutate(H=as.factor(1))
  }


  temp.plot <-
    labelled.input %>% group_by(H,E) %>%
    ggplot(aes(x=name.cov,y=plot.metric)
    ) %>%
    + geom_point(size=point.size) %>%
    + coord_flip()	%>%
    + ylim(lbound,ubound)

  if (diagnostic!=2 & (scope=="all" | (scope=="average" & average.over!="time" & average.over!="distance"))) {
    temp.plot <- temp.plot %>% + facet_grid(rev.exposure~covariate,labeller=label_wrap_gen(width = label.width, multi_line = TRUE))
  } else if (diagnostic==2 & (scope=="all" | (scope=="average" & average.over!="time" & average.over!="distance"))) {
    temp.plot <- temp.plot %>% + facet_grid(rev.covariate~exposure,labeller=label_wrap_gen(width = label.width, multi_line = TRUE))
  } else if (scope=="recent" | (scope=="average" & (average.over=="time" | average.over=="distance"))) {
    temp.plot <- temp.plot %>% + facet_grid(.~comparison,labeller=label_wrap_gen(width = label.width, multi_line = TRUE))
  }

  if (is.null(is.null(axis.title) & metric=="D")) {
    axis.title <- "Mean Difference"
  } else if (is.null(axis.title) & metric=="SMD") {
    axis.title <- "Standardized Mean Difference"
  }

  final.plot <- temp.plot %>%
    + ylab(axis.title) %>%
    + xlab("Covariate") %>%
    + geom_hline(yintercept=c(refline.limit.a,refline.limit.b),linetype="dotted",colour=alpha("black",1/3),size=refline.size) %>%
    + geom_hline(yintercept=0,linetype="solid",colour=alpha("black",.5),size=zeroline.size) %>%
    + scale_x_discrete(limits=rev(sort(unique(input$name.cov)))) %>%
    + themes

  final.plot

}