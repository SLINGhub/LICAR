options(java.parameters = "-Xmx2g" )

library('enviPat') #M calculation for isotopic correction
library('stringr')

##' @title Relative abundance of M+2 of the molecular species with lower molecular mass.
##' @param molec Formula of the lipid without Head Group.
##' @return Relative abundance of M+2.

mCalc <- function (molec)
{
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


##' @title Isotopic correction of the peak from HILIC.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @param lipidGroup The analysis group of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect <- function(inputData, lipidClass, lipidGroup) {
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
  } else { stop("Wrong lipid group! Lipid group should be 'Head Group', 'FA' or 'LCB'.") }
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
  
  if( !(lipidClass %in% c("LPC", "LPCO", "LPCOql", "PC", "SM", "LPE", "PE", "PI", "PG", "PS", "S1P", 
                          "S1Pql", "LPENHG", "PENHG", "PINHG", "PGNHG", "PSNHG") ) ) {
    stop("For lipid group 'Head Group', lipid class must be 'LPC', 'LPCO', 'LPCOql', 'PC', 'SM', 'LPE', 
         'PE', 'PI', 'PG', 'PS', 'S1P', 'S1Pql', 'LPENHG', 'PENHG', 'PINHG', 'PGNHG' or 
         'PSNHG'.") } else {
           
           ##Calculate C, H, O and N, based on the name of lipids
           if(  ( lipidClass %in% c("LPC") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPCO") ) && ( inputData$Product <= 104.3 && inputData$Product >= 103.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 1, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("LPCOql") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 0, constant_O = 2, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PC") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -4, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("SM") ) && ( inputData$Product <= 184.3 && inputData$Product >= 183.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 1 )
           } else if(  ( lipidClass %in% c("LPE") ) && ( ( inputData$Precursor - inputData$Product ) <= 141.2 && ( inputData$Precursor - inputData$Product ) >= 140.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -1, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PE") ) && ( ( inputData$Precursor - inputData$Product ) <= 141.2 && ( inputData$Precursor - inputData$Product ) >= 140.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PI") ) && ( ( inputData$Precursor - inputData$Product ) <= 277.2 && ( inputData$Precursor - inputData$Product ) >= 276.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PG") ) && ( ( inputData$Precursor - inputData$Product ) <= 189.2 && ( inputData$Precursor - inputData$Product ) >= 188.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PS") ) && ( ( inputData$Precursor - inputData$Product ) <= 185.2 && ( inputData$Precursor - inputData$Product ) >= 184.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -3, constant_O = 4, constant_N = 0 )
           } else if(  ( lipidClass %in% c("S1P") ) && ( inputData$Product <= 60.3 && inputData$Product >= 59.9 ) ) {
             isoCorrect_head(inputData, constant_C = 1, constant_H = 1, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("S1Pql") ) && ( inputData$Product <= 113.2 && inputData$Product >= 112.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = 1, constant_O = 1, constant_N = 1 )
           } else if(  ( lipidClass %in% c("LPENHG") ) && ( inputData$Product <= 196.3 && inputData$Product >= 195.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = 0, constant_O = 2, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PENHG") ) && ( inputData$Product <= 196.3 && inputData$Product >= 195.9 ) ) {
             isoCorrect_head(inputData, constant_C = 0, constant_H = -2, constant_O = 3, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PINHG") ) && ( inputData$Product <= 241.2 && inputData$Product >= 240.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PGNHG") ) && ( inputData$Product <= 153.2 && inputData$Product >= 152.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 5, constant_N = 0 )
           } else if(  ( lipidClass %in% c("PSNHG") ) && ( ( inputData$Precursor - inputData$Product ) <= 87.2 && ( inputData$Precursor - inputData$Product ) >= 86.8 ) ) {
             isoCorrect_head(inputData, constant_C = 3, constant_H = -2, constant_O = 8, constant_N = 0 )
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


##' @title Isotopic correction of the peak from HILIC, LCB.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_LCB <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("dhCer", "Cer", "Hex1Cer", "Hex2Cer") ) ) {
    stop("For lipid group 'LCB', lipid class must be 'dhCer', 'Cer', 'Hex1Cer' or 'Hex2Cer'.") 
  } else {
    ##Calculate C, H, O and N, based on the name of lipids
    inputData.annot <- CH_raw(inputData)
    inputData.front <- data.frame(inputData, inputData.annot$C_front, inputData.annot$H_front, check.names = FALSE)
    colnames(inputData.front)[ncol(inputData.front)-1] <- "C_raw"
    colnames(inputData.front)[ncol(inputData.front)] <- "H_raw"
    inputData.back <- data.frame(inputData, inputData.annot$C_back, inputData.annot$H_back, check.names = FALSE)
    colnames(inputData.back)[ncol(inputData.back)-1] <- "C_raw"
    colnames(inputData.back)[ncol(inputData.back)] <- "H_raw"
    
    if( lipidClass %in% "dhCer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 2, constant_O = 1, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 0, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = 2, constant_O = 3, constant_N = 0)$K
    } else if( lipidClass %in% "Hex1Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 6, constant_H = 0, constant_O = 8, constant_N = 0)$K
    } else if( lipidClass %in% "Hex2Cer"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 0, constant_H = 0, constant_O = 0, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 12, constant_H = -2, constant_O = 13, constant_N = 0)$K
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


##' @title Isotopic correction of the peak from HILIC, FA.
##' @param inputData A matrix of lipid peak, with lipids by row and samples by column.
##' @param lipidClass The class of the lipids.
##' @return The isotopic corrected peak.
##' @author Shanshan Ji

isoCorrect_FA <- function(inputData, lipidClass) {
  
  if( !(lipidClass %in% c("PEP", "PCNFA", "PENFA", "PGNFA", "PINFA", "PSNFA") ) ) {
    stop("For lipid group 'FA', lipid class must be 'PEP', 'PCNFA', 'PENFA', 'PGNFA', 'PINFA' or 'PSNFA'.") 
  } else {
    ##Calculate C, H, O and N, based on the name of lipids
    inputData.annot <- CH_raw(inputData)
    inputData.front <- data.frame(inputData, inputData.annot$C_front, inputData.annot$H_front, check.names = FALSE)
    colnames(inputData.front)[ncol(inputData.front)-1] <- "C_raw"
    colnames(inputData.front)[ncol(inputData.front)] <- "H_raw"
    inputData.back <- data.frame(inputData, inputData.annot$C_back, inputData.annot$H_back, check.names = FALSE)
    colnames(inputData.back)[ncol(inputData.back)-1] <- "C_raw"
    colnames(inputData.back)[ncol(inputData.back)] <- "H_raw"
    
    if( lipidClass %in% "PEP"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 2, constant_H = 2, constant_O = 4, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 3, constant_H = -1, constant_O = 3, constant_N = 0)$K
    } else if( lipidClass %in% "PCNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = 0, constant_O = 8, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "PENFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 5, constant_H = 0, constant_O = 6, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "PGNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 6, constant_H = -1, constant_O = 8, constant_N = 0)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "PINFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 9, constant_H = -3, constant_O = 11, constant_N = 0)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
    } else if( lipidClass %in% "PSNFA"  ) {
      inputData$K_front <- CH_K(inputData.front, constant_C = 6, constant_H = -2, constant_O = 8, constant_N = 1)$K
      inputData$K_back <- CH_K(inputData.back, constant_C = 0, constant_H = -1, constant_O = 2, constant_N = 0)$K
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
      if( inputData$DiffPre[j+1]<=2.2 && inputData$DiffPre[j+1]>=1.8 ) {
        if( inputData$DiffPro[j+1]<=0.2 && inputData$DiffPro[j+1]>=-0.2 ) {
          for ( k in (grep("Product", colnames(inputData.isoCorrect))+1):(ncol(inputData.isoCorrect)-2) ) {
            if( inputData.isoCorrect[i, k]>0 ) {
              inputData.isoCorrect[j+1, k] = inputData.isoCorrect[j+1, k] - inputData.isoCorrect[i, k] * inputData.isoCorrect$K_front[i]
            }
          }
        } else if( inputData$DiffPro[j+1]<=2.2 && inputData$DiffPro[j+1]>=1.8 ) {
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