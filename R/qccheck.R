

#' Similarity checker
#' 
#' @description
#' \code{check_similarity} flags respondents whose response patterns
#' are suspiciously similar.
#' 
#' @param Data A dataframe, ideally no more than 3000 columns (otherwise)
#' the algorithm will take a really long time to run
#' 
#' @param Condition1VarInData A character, identifying the first dimension
#' to the cut the similarity test on (optional)
#'
#' @param Condition2VarInData A character, identifying the second dimension
#' to the cut the similarity test on (optional)
#'
#' @param ZThreshold A numeric, the location in the similarity distribution
#' to flag respondents as being suspiciously similar (a z-score)
#'
#' @export

check_similarity<-function(Data,Condition1VarInData=NULL,Condition2VarInData=NULL,ZThreshold) {

	##########
	## 
	
		#Identify the number of conditions
		if (is.null(Condition1VarInData) && is.null(Condition2VarInData)) NumConditions<-0
		if (!is.null(Condition1VarInData) && is.null(Condition2VarInData)) NumConditions<-1
		if (!is.null(Condition1VarInData) && !is.null(Condition2VarInData)) NumConditions<-2
	
		#identify the group counts and the codes
		
		if (NumConditions==1) {
			tempTable<-table(Data[,Condition1VarInData])
			GroupCounts<-matrix(tempTable)
			ConditionCodes1<-as.numeric(names(tempTable))
			}
			
		if (NumConditions==2) {
			tempTable<-table(Data[,Condition1VarInData],Data[,Condition2VarInData])
			GroupCounts<-matrix(tempTable)
			tempConditionCodes1<-as.numeric(rownames(tempTable))
			tempConditionCodes2<-as.numeric(colnames(tempTable))
			
			ConditionCodes1<-rep(tempConditionCodes1,length(tempConditionCodes2))
			ConditionCodes2<-sort(rep(tempConditionCodes2,length(tempConditionCodes1)))
			}
	
	
		if (NumConditions==0) {
			Data$constant<-1
			Condition1VarInData<-"constant"
			NumConditions<-1
			tempTable<-table(Data[,Condition1VarInData])
			GroupCounts<-matrix(tempTable)
			ConditionCodes1<-as.numeric(names(tempTable))
			}
		
		
	##
	##########


	CitySampleCounts=GroupCounts
	
	library(doParallel)
	library(foreach)

	cl <- makeCluster(detectCores())
	registerDoParallel(cl)

	
	
	StoreProblematicIDs<-foreach (cntry = 1:length(CitySampleCounts)) %dopar% {

		# Reinitialize everything
		Data=Data
		NumConditions=NumConditions
		ConditionCodes1=ConditionCodes1
		Condition1VarInData=Condition1VarInData
		if(exists("ConditionCodes2")) ConditionCodes2=ConditionCodes2
		if(exists("Condition2VarInData")) Condition2VarInData=Condition2VarInData
		
		
		
		library(network)
		library(sna)
	
		if (CitySampleCounts[cntry]>1) {
	
			if (NumConditions==2) tempData<-eval(parse(text=paste("Data[Data$",Condition1VarInData,"%in%ConditionCodes1[cntry] & Data$",Condition2VarInData,"%in%ConditionCodes2[cntry],]",sep="")))
			if (NumConditions==1) tempData<-eval(parse(text=paste("Data[Data$",Condition1VarInData,"%in%ConditionCodes1[cntry],]",sep="")))
	


			StoreDistScores<-matrix(0,ncol=CitySampleCounts[cntry],nrow=CitySampleCounts[cntry])
			NumLoop<-((CitySampleCounts[cntry]^2)/2)-(CitySampleCounts[cntry]/2)
				
			
			rowI = 1
			colI = 2
		
			for (i in 1:(NumLoop)) {
			
				temptable<-table(tempData[rowI,]==tempData[colI,])
				StoreDistScores[rowI,colI]<-temptable[2]/sum(temptable)
				StoreDistScores[colI,rowI]<-temptable[2]/sum(temptable)
				ifelse (colI!=CitySampleCounts[cntry],{colI=colI+1},{rowI=rowI+1; colI=rowI+1})
			
				}
	
			# Do network manipulations
			Distances<-sort(matrix(StoreDistScores)[!(matrix(StoreDistScores)%in%0)])
			BetaDistances<-(Distances-mean(Distances))/sd(Distances)
		
			MinValGT2SD<-min(Distances[which(BetaDistances>=ZThreshold)])
	
			AdjacencyMatrix<-matrix(0,ncol=CitySampleCounts[cntry],nrow=CitySampleCounts[cntry])
			AdjacencyMatrix[StoreDistScores>=MinValGT2SD]<-1
		
			#plot(as.network(AdjacencyMatrix),displaylabels=T)
		
			TempProbIDs<-tempData[-isolates(as.network(AdjacencyMatrix)),1]
			if (length(TempProbIDs)>0) {
				#StoreProblematicIDs[[cntry]]<-TempProbIDs
				return(TempProbIDs)
				rm(TempProbIDs)
				}
			
	
			}
		
		}	
		
		return(unlist(lapply(StoreProblematicIDs,as.character)))

}
 
#' Logic checker
#' 
#' @description
#' \code{check_similarity} flags respondents whose response patterns
#' deviate too far from expected (highly correlated) patterns. 
#' 
#' @param Data A dataframe, ideally no more than 3000 columns (otherwise)
#' the algorithm will take a really long time to run
#' 
#' @param CorrelationThreshold A numeric, the correlation threshold
#'
#' @param ZThreshold A numeric, the location in the similarity distribution
#' to flag respondents as being suspiciously similar (a z-score)
#'
#' @export

check_logic<-function(Data,CorrelationThreshold=.5,ZThreshold=2.3) {

	library(igraph)

	# Initialize internal functions
	findClass<-function(Data) {StoreClass<-matrix(,ncol=1,nrow=dim(Data)[2]); for (i in 1:dim(Data)[2]) StoreClass[i]<-class(Data[,i]); return(StoreClass)}
	
	# Test function
	TestPassFail<-function(x,Direction,CurrMeans) {
		if (!any(is.na(x))) {
			if (Direction=="positive") 
				{ # greater than 1.5x above the mean 
				ifelse (x[1]>(CurrMeans[1]+CurrMeans[1]*-.02) && x[2]>(CurrMeans[2]+CurrMeans[2]*-.02), return(1), return(0))
				}
		
			if (Direction=="negative") 
				{
				ifelse ((x[1]>(CurrMeans[1]+CurrMeans[1]*-.02) && x[2]<(CurrMeans[2]+CurrMeans[2]*-.02)) || (x[1]<(CurrMeans[1]+CurrMeans[1]*-.02) && x[2]>(CurrMeans[2]+CurrMeans[2]*-.02)), 
					return(1), return(0))
				}
			}
		if (any(is.na(x))) return(0)
		}
		
	# Dummy reverser
	Reverser<-function(Dummy) {
		tempDummy<-Dummy
		tempDummy[Dummy%in%1]<-0
		tempDummy[Dummy%in%0]<-1
		return(tempDummy)
		}
	
	##################
	## Main routine ##
	##################
		
		# Retrieve varnames
		VarNames<-colnames(Data)
	
		# Retrieve classes for the vars
		CurrClasses<-findClass(Data)
	
		# Select out the varnames that refer to numeric or integer
		TruncVarNames<-VarNames[CurrClasses=="integer" | CurrClasses=="numeric"]
		
		# Identify variables that we don't want to assess their correlation
		ProhibitedWordsForVars<-c("Timers","ID","Serial","SampleProvider","DV","Sniffer","ProgramVersion",
			"DataCollection","COMP","QTA","QCFlags","WaveID","QuotaDaily","CAWI","City","SurveySettings")
	
			RemoveVars<-unlist(lapply(lapply(strsplit(TruncVarNames,"_",""),function(x) which(x%in%ProhibitedWordsForVars)),length))
		
		if (any(RemoveVars>0)) TruncVarNames<-TruncVarNames[-which(RemoveVars>0)]
		

		# Retrieve the correlations among the numeric/integer variables
		StoreCorr<-cor(Data[,TruncVarNames],use="everything")
	
		# Flag the correlation above which we'll find the relevant pairs of variables to test
		#CorrelationThreshold<-quantile(abs(StoreCorr),probs=seq(0, 1, UseQuantile/100),na.rm=T)[8]

		# Identify the pairs of variables that have correlations above the threshold
		StoreCorrFlag<-StoreCorr
		StoreCorrFlag[abs(StoreCorr)>=CorrelationThreshold]<-1
		StoreCorrFlag[abs(StoreCorr)<CorrelationThreshold]<-0
		colnames(StoreCorrFlag)<-TruncVarNames
		rownames(StoreCorrFlag)<-TruncVarNames
		StoreCorrFlag[upper.tri(StoreCorrFlag)]<-NA
		
		# Get pairs of variables which are the "tests"
		Pairs<-get.edgelist(graph.adjacency(StoreCorrFlag))
			Pairs<-Pairs[-which(Pairs[,1]==Pairs[,2]),]
		
		# Apply the test
		StoreRespondentScore<-matrix(0,ncol=1,nrow=dim(Data)[1])
		for (prs in 1:dim(Pairs)[1]) {
			
			# Identify the means of the current pair of variables
			CurrMeans<-apply(Data[,Pairs[prs,]],2,mean,na.rm=T)
			
			# Identify the directionality of the correlation
			ifelse (StoreCorr[Pairs[prs,1],Pairs[prs,2]]<0,CurrDirection<-"negative",CurrDirection<-"positive")
			#print(CurrDirection)
		
			# Apply the test	
				TEST<-apply(Data[,Pairs[prs,]],1,TestPassFail,CurrDirection,CurrMeans)

			# Add the test results to the respondent level, using the CurrDirection
				StoreRespondentScore<-StoreRespondentScore+TEST * abs(StoreCorr[Pairs[prs,1],Pairs[prs,2]])
				
			# Subtract if they fail
				#StoreRespondentScore<-StoreRespondentScore-Reverser(TEST)*20
			}
			

			
		# Flag the people in the bottom part of the distribution as identified by z threshold
			
			StandFlagThreshold<-(StoreRespondentScore-mean(StoreRespondentScore))/sd(StoreRespondentScore)
			
			StandFlagThreshold[which(StandFlagThreshold>2)]<-0 # ignore outliers from the top part of the by setting them to 0
			
			# restandardize
			StandFlagThreshold<-(StandFlagThreshold-mean(StandFlagThreshold))/sd(StandFlagThreshold)
			
			ScoresLessThanZThresh<- StoreRespondentScore[which(StandFlagThreshold< -1*ZThreshold)]
			ifelse (length(ScoresLessThanZThresh)>0,  
				{	
				# Flag and Return the respondent IDs	
				FlagThreshold<-max(ScoresLessThanZThresh)			
				return(Data[StoreRespondentScore<=FlagThreshold,1])
				},
				
				return("NO RESPONDENTS FLAGGED"))
}
