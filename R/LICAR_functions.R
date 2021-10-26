options(java.parameters = "-Xmx2g" )

require('enviPat') #M calculation for isotopic correction
require('stringr')

##' @title Relative abundance of M+2 of the molecular species with lower molecular mass.
##' @param molec Formula of the lipid without Head Group.
##' @return Relative abundance of M+2.

mCalc <- function (molec)
{
  require('enviPat') #M calculation for isotopic correction
  data(isotopes) #Needed in 'enviPat'
  data(resolution_list) #Needed in 'enviPat'
  checked<-check_chemform( isotopes, molec )
  
  centro<-isowrap( isotopes, checked, resmass=resolution_list[[7]], resolution=FALSE, 
                   nknots=4, spar=0.2, threshold=0.1, charge=1, emass=0.00054858, algo=2, ppm=FALSE, 
                   dmz="get", frac=1/4, env="Gaussian", detect="centroid",
                   plotit=FALSE #Plot the grahp if plotit=TRUE
  )
  
  centro <- data.frame(centro)
  K <- centro[3, 2]/100 #M+2
  return(K)
}


mCalc3 <- function (molec)
{
  require('enviPat') #M calculation for isotopic correction
  #data(adducts) 
  data(isotopes) #Needed in 'enviPat'
  data(resolution_list) #Needed in 'enviPat'
  checked<-check_chemform( isotopes, molec )
  
  centro<-isowrap( isotopes, checked, resmass=resolution_list[[7]], resolution=FALSE, 
                   nknots=4, spar=0.2, threshold=0.1, charge=1, emass=0.00054858, algo=2, ppm=FALSE, 
                   dmz="get", frac=1/4, env="Gaussian", detect="centroid",
                   plotit=FALSE #Plot the grahp if plotit=TRUE
  )
  
  centro <- data.frame(centro)
  M3 <- centro[4, 2]/100 #M+3
  return(M3)
}

##################################################################################Main functions
##' @title Isotopic correction of the peak from HILIC.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @param lipidGroup The analysis group of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect <- function(inputData, lipidClass, lipidGroup) {
  
  colnames(inputData)[1] <- "Precursor"
  colnames(inputData)[2] <- "Product"
  
  if( length(inputData$Precursor) == 0 ) stop("Precursor is missing!")
  if( length(inputData$Product) == 0 ) stop("Product is missing!")
  
  inputData <- inputData[order(inputData$Precursor, inputData$Product), ] #Sort the inputData
  inputData[is.na(inputData)] <- 0 ##Replace NA by 0
  
  #Check the lipidGroup
  if(lipidGroup %in% "Head Group") {
    isoCorrect_headGroup(inputData, lipidClass)
  } else if(lipidGroup %in% "FA") {
    isoCorrect_FA(inputData, lipidClass)
  } else if(lipidGroup %in% "LCB") {
    isoCorrect_LCB(inputData, lipidClass) 
  } else if(lipidGroup %in% "Neutral") {
    isoCorrect_Neutral(inputData, lipidClass)
  } else if(lipidGroup %in% "RPLC") {
    isoCorrect_RPLC(inputData, lipidClass)
  } else { stop("Wrong lipid group! Lipid group should be 'Head Group', 'FA', 'LCB', 'Neutral' or 'RPLC'.") }
}


##' @title Isotopic correction of the peak from HILIC, Pos (head group) or Neg (head group).
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param constant_C Constant for C.
##' @param constant_H Constant for H.
##' @param constant_O Constant for O.
##' @param constant_N Constant for N, 0 or 1.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_head <- function(inputData, constant_C, constant_H, constant_O, constant_N) {
  ##Get the position of ":"
  position <- regexpr(pattern =':', rownames(inputData)) #Returns position of 1st match in a string
  position <- matrix(position)
  
  ##Calculate C, H, O and N, based on the name of lipids
  inputData$C_raw <- as.numeric ( substring(rownames(inputData), position - 2, position - 1) )
  inputData$H_raw <- as.numeric ( substring(rownames(inputData), position + 1, position + 2) ) 
  inputData$C <- inputData$C_raw + constant_C
  inputData$H <- inputData$C * 2 + constant_H - 2 * inputData$H_raw
  
  if(constant_N == 0 && constant_O == 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, sep="")
  } else if (constant_N == 0 && constant_O > 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "O", constant_O, sep="")
  } else if (constant_N > 0 && constant_O == 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "N", constant_N, sep="")
  } else if (constant_N > 0 && constant_O > 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "N", constant_N, "O", constant_O, sep="")
  }
  
  ##Get the relative abundance
  inputData$K <- 0
  
  for ( i in 1:nrow(inputData) )
  {
    inputData$K[i] <- mCalc(inputData$Formula[i])
  } 
  
  ##Check the difference of precursor
  inputData$Diff <- 0
  
  for ( i in 2:nrow(inputData) )
  {
    inputData$Diff[i] <- inputData$Precursor[i] - inputData$Precursor[i-1]
  }
  
  ##Isotopic correction, Intensity'[i] = Intensity[i] - K[i-1] * Intensity[i-1]
  inputData.isoCorrect <- inputData
  #dim(inputData.isoCorrect)
  ########Recalculate the peak
  for ( i in 2:nrow(inputData.isoCorrect))
  {
    if( inputData.isoCorrect$Diff[i] <= 2.2 && inputData.isoCorrect$Diff[i] >= 1.8 )  {
      for ( j in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect) - 7) ) {
        if ( inputData.isoCorrect[i-1, j] > 0 ) {
          inputData.isoCorrect[i, j] = inputData.isoCorrect[i, j] - inputData.isoCorrect[i-1, j] *  inputData.isoCorrect$K[i-1] 
        }
      }
    }
  }
  inputData.isoCorrect[inputData.isoCorrect<0] <- 0 #Replace the negtive numbers by 0
  inputData.isoCorrect <- inputData.isoCorrect[, c(1:(ncol(inputData.isoCorrect)-7))]
  return(inputData.isoCorrect) 
} 


##' @title Isotopic correction of the peak from HILIC, head group.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_headGroup <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("AcylCarnitine", "LPC", "LPCql", "LPCO", "LPCOql", "LPE", "LPENHG", "PC", "PCO", "PCP", "PE", "PENHG", 
                          "PG", "PGNHG", "PI", "PINHG", "PS", "PSNHG", "S1P", "S1Pql", "SM") ) ) {
    stop("For lipid group 'Head Group', lipid class must be 'AcylCarnitine', 'LPC', 'LPCql', 'LPCO', 'LPCOql', 'LPE', 'LPENHG', 'PC', 'PCO', 'PCP', 'PE', 'PENHG', 'PG', 
         'PGNHG', 'PI', 'PINHG', 'PS', 'PSNHG', 'S1P', 'S1Pql' or 'SM'.") } else {
           
           ##Calculate C, H, O and N, based on the name of lipids
           if(  ( lipidClass %in% c("AcylCarnitine") ) && ( inputData$Product <= 85.2 && inputData$Product >= 84.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 3, constant_O = 2, constant_N = 1 )
           } else if(  ( lipidClass %in% c("LPC") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPCql") ) && ( inputData$Product <= 104.3 && inputData$Product >= 103.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -1, constant_O = 6, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPCO") ) && ( inputData$Product <= 104.3 && inputData$Product >= 103.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 1, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPCOql") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 0, constant_O = 2, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPE") ) && ( ( inputData$Precursor - inputData$Product ) <= 141.2 && ( inputData$Precursor - inputData$Product ) >= 140.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -1, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPENHG") ) && ( inputData$Product <= 196.3 && inputData$Product >= 195.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = 0, constant_O = 2, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PC") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -4, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PCO") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PCP") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -4, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PE") ) && ( ( inputData$Precursor - inputData$Product ) <= 141.2 && ( inputData$Precursor - inputData$Product ) >= 140.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PENHG") ) && ( inputData$Product <= 196.3 && inputData$Product >= 195.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = -2, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PG") ) && ( ( inputData$Precursor - inputData$Product ) <= 189.2 && ( inputData$Precursor - inputData$Product ) >= 188.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PGNHG") ) && ( inputData$Product <= 153.2 && inputData$Product >= 152.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PI") ) && ( ( inputData$Precursor - inputData$Product ) <= 277.2 && ( inputData$Precursor - inputData$Product ) >= 276.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PINHG") ) && ( inputData$Product <= 241.2 && inputData$Product >= 240.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PS") ) && ( ( inputData$Precursor - inputData$Product ) <= 185.2 && ( inputData$Precursor - inputData$Product ) >= 184.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PSNHG") ) && ( ( inputData$Precursor - inputData$Product ) <= 87.2 && ( inputData$Precursor - inputData$Product ) >= 86.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 8, constant_N = 0 )
           } else if(  ( lipidClass %in% c("S1P") ) && ( inputData$Product <= 60.3 && inputData$Product >= 59.9 ) ) {
             isoCorrect_head(inputData, constant_C = 1, constant_H = 1, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("S1Pql") ) && ( inputData$Product <= 113.2 && inputData$Product >= 112.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 1, constant_O = 1, constant_N = 1 )
           } else if(  ( lipidClass %in% c("SM") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 1 )
           }
         }  
} 


##' @title Calculate C and H based on the first and second ":".
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @return The number for C and H calculation.
##' @author Shanshan Ji

CH_raw <- function(inputData) {
  ##Get the C_front and H_front by the first ":"
  position = regexpr(pattern =':', rownames(inputData)) #Returns position of 1st match in a string
  position = matrix(position)
  inputData$C_front = as.numeric ( substring(rownames(inputData), position - 2, position - 1) )
  inputData$H_front = as.numeric ( substring(rownames(inputData), position + 1, position + 1) ) 
  
  #Get the formula2 for FA
  k=strsplit(rownames(inputData), ":\\s*(?=[^:]+$)", perl=TRUE)
  l=do.call(mapply, c(cbind, k))
  inputData$C_back = as.numeric ( str_sub(l[, 1], -2, -1) )
  inputData$H_back = as.numeric ( str_sub(l[, 2], 1, 1) ) 
  
  return(inputData[, c((ncol(inputData)-3):ncol(inputData))]) 
}


##' @title Calculate the relative abundance for the isotopic correction.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param constant_C Constant for C.
##' @param constant_H Constant for H.
##' @param constant_O Constant for O.
##' @param constant_N Constant for N, 0 or 1.
##' @return The relative abundance.
##' @author Shanshan Ji

CH_K <- function(inputData, constant_C, constant_H, constant_O, constant_N) {
  inputData$C <- inputData$C_raw + constant_C
  inputData$H <- inputData$C * 2 + constant_H - 2 * inputData$H_raw
  
  if(constant_N == 0 && constant_O == 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, sep="")
  } else if (constant_N == 0 && constant_O > 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "O", constant_O, sep="")
  } else if (constant_N > 0 && constant_O == 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "N", constant_N, sep="")
  } else if (constant_N > 0 && constant_O > 0) {
    inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "N", constant_N, "O", constant_O, sep="")
  }
  
  ##Get the relative abundance
  inputData$K <- 0
  
  for ( i in 1:nrow(inputData) )
  {
    inputData$K[i] <- mCalc(inputData$Formula[i])
  } 
  return(inputData) 
}  



##' @title Isotopic correction of the peak from HILIC, FA.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_FA <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("CLNFA", "CLNFA_2", "LPCNFA", "LPCNFA_2", "LPENFA", "LPINFA", "LPGNFA", 
                          "PCNFA", "PCONFA", "PCPNFA", "PCNFA_2", "PCONFA_2", "PCPNFA_2", 
                          "PENFA", "PEONFA", "PEP", "PENFA", "PGNFA", "PINFA", "PSNFA") ) ) {
    stop("For lipid group 'FA', lipid class must be 'CLNFA', 'CLNFA_2', 'LPCNFA', 'LPCNFA_2', 
         'LPENFA', 'LPINFA', 'LPGNFA', 'PCNFA', 'PCONFA', 'PCPNFA', 'PCNFA_2', 'PCONFA_2', 'PCPNFA_2', 
         'PENFA', 'PEONFA', 'PEP', 'PENFA', 'PGNFA', 'PINFA' or 'PSNFA'.") 
  } else {
    ##Calculate C, H, O and N, based on the name of lipids
    inputData.annot <- CH_raw(inputData)
    inputData.front <- data.frame(inputData, inputData.annot$C_front, inputData.annot$H_front, check.names = FALSE)
    colnames(inputData.front)[ncol(inputData.front)-1] <- "C_raw"
    colnames(inputData.front)[ncol(inputData.front)] <- "H_raw"
    inputData.back <- data.frame(inputData, inputData.annot$C_back, inputData.annot$H_back, check.names = FALSE)
    colnames(inputData.back)[ncol(inputData.back)-1] <- "C_raw"
    colnames(inputData.back)[ncol(inputData.back)] <- "H_raw"
    
    if(  ( lipidClass %in% c("LPCNFA") ) && ( ( inputData$Precursor - inputData$Product ) <= 285.3 && ( inputData$Precursor - inputData$Product ) >= 284.9 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0 )
      return(inputData.isoCorrect)
    } else if(  ( lipidClass %in% c("LPCNFA_2") ) && ( ( inputData$Precursor - inputData$Product ) <= 299.3 && ( inputData$Precursor - inputData$Product ) >= 298.9 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0 )
      return(inputData.isoCorrect)
    } else if(  ( lipidClass %in% c("LPENFA") ) && ( ( inputData$Precursor - inputData$Product ) <= 197.3 && ( inputData$Precursor - inputData$Product ) >= 196.9 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0 )
      return(inputData.isoCorrect)
    } else if(  ( lipidClass %in% c("LPINFA") ) && ( ( inputData$Precursor - inputData$Product ) <= 316.3 && ( inputData$Precursor - inputData$Product ) >= 315.9 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0 )
      return(inputData.isoCorrect)
    } else if(  ( lipidClass %in% c("LPGNFA") ) && ( ( inputData$Precursor - inputData$Product ) <= 228.3 && ( inputData$Precursor - inputData$Product ) >= 227.9 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0 )
      return(inputData.isoCorrect)
    } else if( lipidClass %in% "PCNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = 0, constant_O = 8, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PCONFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = 2, constant_O = 7, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PCPNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = 0, constant_O = 7, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PCNFA_2"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 10, constant_H = 0, constant_O = 8, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PCONFA_2"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 10, constant_H = 2, constant_O = 7, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PCPNFA_2"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 10, constant_H = 0, constant_O = 7, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PENFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 5, constant_H = 0, constant_O = 6, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PEONFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 5, constant_H = 2, constant_O = 5, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PEP"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 2, constant_H = 2, constant_O = 4, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 3, constant_H = -1, constant_O = 3, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PGNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 6, constant_H = -1, constant_O = 8, constant_N = 0)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PINFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = -3, constant_O = 11, constant_N = 0)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "PSNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 6, constant_H = -2, constant_O = 8, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "CLNFA"  ) { #"CLNFA" is special, so K_front is calculated separately
      inputData$C <- inputData.annot$C_front - inputData.annot$C_back + 9
      inputData$H <- inputData$C * 2 - 4 - 2 * ( inputData.annot$H_front - inputData.annot$H_back )
      
      inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "O15", sep="")
      
      ##Get the relative abundance
      inputData$K_front <- 0
      
      for ( i in 1:nrow(inputData) )
      {
        inputData$K_front[i] <- mCalc(inputData$Formula[i])
      } 
      
      inputData <- inputData[, !(colnames(inputData) %in% c("C", "H", "Formula"))]
      
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      
      DiffPre_front <- 2
      DiffPro_front <- 0
      DiffPro_back <- 2
    } else if( lipidClass %in% "CLNFA_2"  ) { #"CLNFA" is special, so K_front is calculated separately
      inputData$C <- inputData.annot$C_front - inputData.annot$C_back + 9
      inputData$H <- inputData$C * 2 - 4 - 2 * ( inputData.annot$H_front - inputData.annot$H_back )
      
      inputData$Formula <- paste("C", inputData$C, "H", inputData$H, "O15", sep="")
      
      ##Get the relative abundance
      inputData$K_front <- 0
      
      for ( i in 1:nrow(inputData) )
      {
        inputData$K_front[i] <- mCalc(inputData$Formula[i])
      } 
      
      inputData <- inputData[, !(colnames(inputData) %in% c("C", "H", "Formula"))]
      
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
      
      DiffPre_front <- 1
      DiffPro_front <- 0
      DiffPro_back <- 2
    } 
  }
  
  
  ##Isotopic correction, Intensity'[i] = Intensity[i] - Intensity[i-1] * K[i-1]
  inputData.isoCorrect <- inputData
  ##Check the difference of precursor
  #The inputData are sorted by Precursor first then Product already
  for ( i in 1:( nrow(inputData)-1 ) ) {
    inputData$DiffPre <- 0
    inputData$DiffPro <- 0
    for ( j in i:( nrow(inputData)-1) ) {
      inputData$DiffPre[j+1] <- inputData$Precursor[j+1] - inputData$Precursor[i]
      inputData$DiffPro[j+1] <- inputData$Product[j+1] - inputData$Product[i]
      if( inputData$DiffPre[j+1]<=DiffPre_front + 0.2 && inputData$DiffPre[j+1]>= DiffPre_front - 0.2 ) {
        if( inputData$DiffPro[j+1]<= DiffPro_front + 0.2 && inputData$DiffPro[j+1]>= DiffPro_front - 0.2 ) {
          for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
            if( inputData.isoCorrect[i, k]>0 ) {
              inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_front[i]
            }
          }
        } else if( inputData$DiffPro[j+1]<= DiffPro_back + 0.2 && inputData$DiffPro[j+1]>= DiffPro_back - 0.2 ) {
          for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
            if( inputData.isoCorrect[i, k]>0 ) {
              inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_back[i]
            }
          }
        }
      }
    }
  }
  inputData.isoCorrect <- inputData.isoCorrect[, 1:(ncol(inputData.isoCorrect)-2)]
  inputData.isoCorrect[ inputData.isoCorrect<0 ] <- 0
  return(inputData.isoCorrect) 
}




##' @title Isotopic correction of the peak from HILIC, LCB.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_LCB <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("Cer", "deoxyCer", "dhCer", "dhCer_2", "GM3", "Hex1Cer", "Hex2Cer", "Hex3Cer", "Hex1Sph") ) ) {
    stop("For lipid group 'LCB', lipid class must be 'Cer', 'deoxyCer', 'dhCer', 'dhCer_2', 'GM3', 'Hex1Cer', 'Hex2Cer', 'Hex3Cer' or 'Hex1Sph'.") 
  } else {
    ##Calculate C, H, O and N, based on the name of lipids
    inputData.annot <- CH_raw(inputData)
    inputData.front <- data.frame(inputData, inputData.annot$C_front, inputData.annot$H_front, check.names = FALSE)
    colnames(inputData.front)[ncol(inputData.front)-1] <- "C_raw"
    colnames(inputData.front)[ncol(inputData.front)] <- "H_raw"
    inputData.back <- data.frame(inputData, inputData.annot$C_back, inputData.annot$H_back, check.names = FALSE)
    colnames(inputData.back)[ncol(inputData.back)-1] <- "C_raw"
    colnames(inputData.back)[ncol(inputData.back)] <- "H_raw"
    
    if(  ( lipidClass %in% c("Hex1Sph") ) &&( ( inputData$Precursor - inputData$Product ) <= 180.2 && ( inputData$Precursor - inputData$Product ) >= 179.8 ) ) {
      inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = 2, constant_O = 1, constant_N = 1 )
      return(inputData.isoCorrect)
    } else if( lipidClass %in% "Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 2, constant_O = 3, constant_N = 0)$K
    } else if( lipidClass %in% "deoxyCer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 2, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 0, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "dhCer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 2, constant_O = 1, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 0, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "dhCer_2"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 2, constant_O = 3, constant_N = 0)$K
    } else if( lipidClass %in% "GM3"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 23, constant_H = -7, constant_O = 21, constant_N = 1)$K
    } else if( lipidClass %in% "Hex1Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 6, constant_H = 0, constant_O = 8, constant_N = 0)$K
    } else if( lipidClass %in% "Hex2Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 12, constant_H = -2, constant_O = 13, constant_N = 0)$K
    } else if( lipidClass %in% "Hex3Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 18, constant_H = -4, constant_O = 18, constant_N = 0)$K
    } 
  }
  
  ##Isotopic correction, Intensity'[i] = Intensity[i] - Intensity[i-1] * K[i-1]
  inputData.isoCorrect <- inputData
  ##Check the difference of precursor
  for ( i in 1:( nrow(inputData)-1 ) ) {
    inputData$DiffPre <- 0
    inputData$DiffPro <- 0
    for ( j in i:( nrow(inputData)-1) ) {
      inputData$DiffPre[j+1] <- inputData$Precursor[j+1] - inputData$Precursor[i]
      inputData$DiffPro[j+1] <- inputData$Product[j+1] - inputData$Product[i]
      if( inputData$DiffPre[j+1]<=2.2 && inputData$DiffPre[j+1]>=1.8 ) {
        if( inputData$DiffPro[j+1]<=0.2 && inputData$DiffPro[j+1]>=-0.2 ) {
          for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
            if( inputData.isoCorrect[i, k]>0 ) {
              inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_back[i]
            }
          }
        } else if( inputData$DiffPro[j+1]<=2.2 && inputData$DiffPro[j+1]>=1.8 ) {
          for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
            if( inputData.isoCorrect[i, k]>0 ) {
              inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_front[i]
            }
          }
        }
      }
    }
  }
  inputData.isoCorrect <- inputData.isoCorrect[, 1:(ncol(inputData.isoCorrect)-2)]
  inputData.isoCorrect[ inputData.isoCorrect<0 ] <- 0
  return(inputData.isoCorrect) 
}


##' @title Isotopic correction of the peak from HILIC, head group.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_Neutral <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("CE", "DG", "TG") ) ) {
    stop("For lipid group 'Neutral', lipid class must be 'CE', 'DG' or 'TG'.") } else {
      
      if( ( lipidClass %in% c("CE") ) && ( inputData$Product <= 369.5 && inputData$Product >= 369.1 ) ) {
        inputData.isoCorrect <- isoCorrect_head(inputData, constant_C = 0, constant_H = 3, constant_O = 2, constant_N = 1 )
      }  else if( lipidClass %in% c("DG", "TG") ) {
        
        ##Calculate C, H, O and N, based on the name of lipids
        inputData.annot <- CH_raw(inputData)
        
        #The inputData.front is inputData.mix in fact, espcially for DG,TG of Neutral
        inputData.front <- data.frame(inputData, inputData.annot$C_front - inputData.annot$C_back, 
                                      inputData.annot$H_front - inputData.annot$H_back, 
                                      check.names = FALSE)
        colnames(inputData.front)[ncol(inputData.front)-1] <- "C_raw"
        colnames(inputData.front)[ncol(inputData.front)] <- "H_raw"
        
        inputData.back <- data.frame(inputData, inputData.annot$C_back, inputData.annot$H_back, check.names = FALSE)
        colnames(inputData.back)[ncol(inputData.back)-1] <- "C_raw"
        colnames(inputData.back)[ncol(inputData.back)] <- "H_raw"
        
        
        if( lipidClass %in% "DG"  ) {
          inputData$K_front <- CH_K(inputData.front, constant_C = 3, constant_H = -1, constant_O = 3, constant_N = 0)$K
          inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 3, constant_O = 2, constant_N = 1)$K
        } else if( lipidClass %in% "TG"  ) {
          inputData$K_front <- CH_K(inputData.front, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0)$K
          inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 3, constant_O = 2, constant_N = 1)$K
        } 
        
        ##Isotopic correction, Intensity'[i] = Intensity[i] - Intensity[i-1] * K[i-1]
        inputData.isoCorrect <- inputData
        ##Check the difference of precursor
        for ( i in 1:( nrow(inputData)-1 ) ) {
          inputData$DiffPre <- 0
          inputData$DiffPro <- 0
          for ( j in i:( nrow(inputData)-1) ) {
            inputData$DiffPre[j+1] <- inputData$Precursor[j+1] - inputData$Precursor[i]
            inputData$DiffPro[j+1] <- inputData$Product[j+1] - inputData$Product[i]
            if( inputData$DiffPre[j+1]<=2.2 && inputData$DiffPre[j+1]>=1.8 ) {
              if( inputData$DiffPro[j+1]<=0.2 && inputData$DiffPro[j+1]>=-0.2 ) {
                for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
                  if( inputData.isoCorrect[i, k]>0 ) {
                    inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_back[i]
                  }
                }
              } else if( inputData$DiffPro[j+1]<=2.2 && inputData$DiffPro[j+1]>=1.8 ) {
                for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
                  if( inputData.isoCorrect[i, k]>0 ) {
                    inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_front[i]
                  }
                }
              }
            }
          }
        }
        inputData.isoCorrect <- inputData.isoCorrect[, 1:(ncol(inputData.isoCorrect)-2)]
        inputData.isoCorrect[ inputData.isoCorrect<0 ] <- 0
        
      } 
      return(inputData.isoCorrect) 
    }       
}



#########################################################################################################################################

##' @title Isotopic correction of the peak: PC by SM and PC-O by PC-P.
##' @param data A matrix of isotopic corrected lipid peak, lipids by row and samples by column.
##' @param class The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_RPLC <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("PC", "PCO") ) ) {
    stop("For lipid group 'Head Group', lipid class must be 'PC' or 'PCO'.") } else {
      
      ##Replace NA by 0
      inputData[is.na(inputData)] <- 0
      
      ##Sort by Precursor and	Product
      inputData <- inputData[order(inputData$Precursor, inputData$Product), ]
      
      ##Get the position of ":"
      position <- regexpr(pattern =':', rownames(inputData)) #Returns position of 1st match in a string
      position <- matrix(position)
      
      ##Calculate Cmatch and Hmatch
      inputData$Cmatch <- as.numeric ( substring(rownames(inputData), position - 2, position - 1) )
      inputData$Hmatch <- as.numeric ( substring(rownames(inputData), position + 1, position + 1) ) 
      
      ##Get class
      inputData$Class <- gsub(" .*$", "", rownames(inputData))
      inputData.PC <- inputData[inputData$Class %in% "PC", ]
      inputData.PCO <- inputData[inputData$Class %in% "PC-O", ]
      inputData.PCP <- inputData[inputData$Class %in% "PC-P", ]
      inputData.SM <- inputData[inputData$Class %in% "SM", ]
      
      ###################################PC corrected by SM
      if (!is.na(inputData.SM[1, 1])) {       #if (!(dim(inputData.SM)[1] == 0)) {
        
        ##Calculate C and H 
        inputData.SM$C <- inputData.SM$Cmatch
        inputData.SM$H <- inputData.SM$C * 2 - 1 - 2 * inputData.SM$Hmatch
        
        ##Get the formula
        inputData.SM$Formula <- paste("C", inputData.SM$C, "H", inputData.SM$H, "NO2", sep="")
        
        inputData.SM$M <- 0
        
        for ( i in 1:nrow(inputData.SM) )
        {
          inputData.SM$M[i] <- mCalc3(inputData.SM$Formula[i])
        } 
        
        #Correction
        inputData.PC.new <- inputData.PC
        #colnames(inputData.PC)
        for ( j in 3:(ncol(inputData.PC)-3) ) 
        {
          for ( i in 1:nrow(inputData.PC) )
          {
            for ( k in 1:nrow(inputData.SM) )
            {
              if ( (inputData.SM$Cmatch[k] - inputData.PC$Cmatch[i] == 4) && (inputData.SM$Hmatch[k] - inputData.PC$Hmatch[i] == 1) )
              {
                inputData.PC.new[i, j] <- inputData.PC[i, j] - inputData.SM[k, j] *  inputData.SM$M[k]
              }
            }
          }
        }
        
        ##Replace the number<=0 by the 0.1
        inputData.PC.new[, 3:(ncol(inputData.PC.new)-3)][inputData.PC.new[, 3:(ncol(inputData.PC.new)-3)] <= 0] <- 0.1
        
        inputData.PC_SM.new <- rbind(inputData.PC.new[, 1:ncol(inputData)], inputData.SM[, 1:ncol(inputData)])
        
      }
      
      ###################################PCO corrected by PCP
      ##Calculate C and H 
      if (!is.na(inputData.PCO[1, 1])) {         #if (!(dim(inputData.PCO)[1] == 0)) {
        
        inputData.PCP$C <- inputData.PCP$Cmatch + 3
        inputData.PCP$H <- inputData.PCP$C * 2 - 4 - 2 * inputData.PCP$Hmatch
        
        ##Get the formula
        inputData.PCP$Formula <- paste("C", inputData.PCP$C, "H", inputData.PCP$H, "O3", sep="")
        
        inputData.PCP$M <- 0
        
        for ( i in 1:nrow(inputData.PCP) )
        {
          inputData.PCP$M[i] <- mCalc(inputData.PCP$Formula[i])
        } 
        
        #Correction
        inputData.PCO.new <- inputData.PCO
        #colnames(inputData.PCO)
        for ( j in 3:(ncol(inputData.PCO)-3) ) 
        {
          for ( i in 1:nrow(inputData.PCO) )
          {
            for ( k in 1:nrow(inputData.PCP) )
            {
              if ( (inputData.PCP$Cmatch[k] == inputData.PCO$Cmatch[i] ) && (inputData.PCP$Hmatch[k] == inputData.PCO$Hmatch[i] ) )
              {
                inputData.PCO.new[i, j] <- inputData.PCO[i, j] - inputData.PCP[k, j] *  inputData.PCP$M[k]
              }
            }
          }
        }
        
        ##Replace the number<=0 by 0.1
        inputData.PCO.new[, 3:(ncol(inputData.PCO.new)-3)][inputData.PCO.new[, 3:(ncol(inputData.PCO.new)-3)] <= 0] <- 0.1
        
        inputData.PCO_PCP.new <- rbind(inputData.PCO.new[, 1:ncol(inputData)], inputData.PCP[, 1:ncol(inputData)])
        
      }
      
      ###Combine the inputData
      if (is.na(inputData.PCO[1, 1])) {
        inputData.new <- inputData.PC_SM.new 
      } else if (is.na(inputData.PC[1, 1])) {
        inputData.new <- inputData.PCO_PCP.new 
      } else inputData.new <- rbind(inputData.PC_SM.new, inputData.PCO_PCP.new)
      
      inputData.new <- inputData.new[, -((ncol(inputData.new)-2):ncol(inputData.new))]
      inputData.new <- inputData.new[order(inputData.new$Precursor, inputData.new$Product), ]
      return(inputData.new) 
    }
}



