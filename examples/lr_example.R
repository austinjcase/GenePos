###########################################################################
#	RGenePos															#######
#	Example Wheat leaf rust from									#######
#	Based on paper Wamishe et al. 2004								#######
#	Version update 21mar16											#######
###########################################################################

## set working directory

setwd("/Users/case0197/Documents/GenePos/example data")

##########################################################################
# load example data from Wamishe et al. 2004 						######
##########################################################################


## this is the data file for unknowns to postulate
lines_IT<-read.csv("lines_IT.csv", head=T, row.names=1)

## this is the data file for single gene differentials
isolines_IT<-read.csv("isolines_IT.csv", head=T, row.names=1)

##########################################################################
# convert raw IT calls                                                 ###
# using the SrIT.Convert functoin to convert IT calls to high "H" and  ###
# low "L" calls basedon the the cutoff value "33" is minimum "L" value ###
# in this example													   ###
##########################################################################

## convert the unknowns
lines_HL<-SrIT.Convert(lines_IT, "33")

## convert the differentials
isolines_HL<-SrIT.Convert(isolines_IT, "33")

##########################################################################
# postulate genes                                                      ###
# using the GenePos function posulate which resistance genes the       ###
# the unknowns might carry. Where "1" is the first unknown in the data ###
# frame of unknowns. Alternatively one may call the unknown by name if ###
# line.name=T is used.
##########################################################################


## postulate the first unknown
GenePos(1, lines_HL, isolines_HL, isolines_IT, lines_IT )

## postulate the second unknown
GenePos(2, lines_HL, isolines_HL, isolines_IT, lines_IT )

## the name of the unknow may also be used
GenePos("MADISON", lines_HL, isolines_HL, isolines_IT, lines_IT, line.name=TRUE)


