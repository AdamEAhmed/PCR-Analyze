#Clear
rm(list = ls())

#Import and combine key file and raw file

##SELECT RAW DATA FILE##
rawCT <- read.table(choose.files(caption = "SELECT RAW DATA FILE"), header = TRUE, sep="\t")

##SELECT KEY FILE##
key <- read.table(choose.files(caption = "SELECT KEY FILE"), header = TRUE, sep="\t")

##SELECT GENE FILE##
gene <- read.table(choose.files(caption = "SELECT GENE FILE"), header = TRUE, sep="\t")

#Merge sets
comp <- merge(rawCT, key, by = 0, sort = FALSE)
comp$Row.names <- NULL
comp$well.y <- NULL

###Set to 1 to omit outliers or set to 0 to include them. (default 1)###
omit <- 1

#Add Delta CT column
comp$Dct <- rep(0,nrow(comp))

#Calculate delta ct
for (z in 1:nrow(gene)) {
  m <- 1
  n <- 2
  for (a in 1:nrow(comp)) {
    if(is.na(comp$gene[m])||is.na(comp$gene[n])){
      m <- m + 1
      n <- n + 1
      next
    }
    #Find control gene of intrest index
    while(comp$gene[m] != gene[z,1]){
      m <- m + 1
      if(m > nrow(comp)){
        break
      }
    }  
    while(comp$gene[n] != gene[z,2]){
      n <- n + 1
      if(n > nrow(comp)){
        break
      }
    }
    if(m > nrow(comp) || n > nrow(comp)){
      break
    }
    testC <- m + 1
    test <- n + 1
    intial <- m
    while (comp$gene[n] == comp$gene[test]) {
      if (comp$gene[m] != comp$gene[testC]){
        comp$Dct[m] <- 0
        comp$Dct[n] <- comp$Ct_value[m] - comp$Ct_value[n]
        m <- intial
        n <- n + 1
        test <- test + 1
        testC <- m + 1
        next
      }
      if (is.na(comp$Ct_value[m]) || is.na(comp$Ct_value[n])){
        comp$Dct[n] <- NA
        n <- n + 1
        m <- m + 1
        test <- test + 1
        testC <- testC + 1
        if (test > nrow(comp)){
          break
        }
        next
      }
      comp$Dct[m] <- 0
      comp$Dct[n] <- comp$Ct_value[m] - comp$Ct_value[n]
      m <- m + 1
      n <- n + 1
      test <- test + 1
      testC <- testC + 1
      if(test > nrow(comp) || n > nrow(comp)){
        break
      }
    }
    comp$Dct[m] <- 0
    comp$Dct[n] <- comp$Ct_value[m] - comp$Ct_value[n]
    m <- m + 1
    n <- n + 1
    if (n > nrow(comp)){
      break
    }
  }
}

#Fixing technical replicates
x <- 1
y <- 2
for (a in 1:nrow(comp)) {
  if(is.na(comp$technical_replicate[x])||is.na(comp$technical_replicate[y])){
    x <- x + 1
    y <- y + 1
    next
  }
  s <- y + 1
  if(is.na(comp$Dct[x]) && comp$technical_replicate[x]==1){
    if(is.na(comp$Dct[y]) && comp$technical_replicate[y]==2){
      comp$technical_replicate[s] <- 1
    }
    comp$technical_replicate[y] <- 1
  }
  else if(comp$technical_replicate[x]==3 && comp$technical_replicate[y]==2){
    comp$technical_replicate[y] <- 1
  }
  x <- x + 1
  y <- y + 1
  if(x > nrow(comp) || y > nrow(comp)){
    break
  }
}

result <- subset(comp, (Ct_value>0)&(Dct>=0|Dct<=0), select = c(1:7))
resultNA <- subset(comp, is.na(Ct_value), select = c(1:6))
result$mean <- NA

#Function to omit outliers
find.out <- function(vec, col){

  count <- 0
  if(length(vec)==3){
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1 &&
       (abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
      vec <- vec[-1]
    }
    else if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1 &&
            (abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
      vec <- vec[-2]
    }
    else if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1 &&
            (abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
      vec <- vec[-3]
    }
  }
  
  if(length(vec)==4){
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(
        vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if (count >= 2){
            vec <- vec[-1]
          }
        }
      }
    }
    if (count >= 2){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if (count >= 2){
            vec <- vec[-2]
          }
        }
      }
    }
    if (count >= 2){
    return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if (count >= 2){
            vec <- vec[-3]
          }
        }
      }
    }
    if (count >= 2){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
          count <- count + 1
          if (count >= 2){
            vec <- vec[-4]
          }
        }
      }
    }
  }
  
  if(length(vec)==5){
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[1])-abs(vec[5]))/((abs(vec[1])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if (count >= 3){
              vec <- vec[-1]
            }
          }
        }
      }
    }
    if (count >= 3){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[2])-abs(vec[5]))/((abs(vec[2])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if (count >= 3){
              vec <- vec[-2]
            }
          }
        }
      }
    }
    if (count >= 3){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[3])-abs(vec[5]))/((abs(vec[3])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if (count >= 3){
              vec <- vec[-3]
            }
          }
        }
      }
    }
    if (count >= 3){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[4])-abs(vec[5]))/((abs(vec[4])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if (count >= 3){
              vec <- vec[-4]
            }
          }
        }
      }
    }
    if (count >= 3){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[5])-abs(vec[2]))/((abs(vec[5])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[5])-abs(vec[3]))/((abs(vec[5])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[5])-abs(vec[4]))/((abs(vec[5])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[5])-abs(vec[1]))/((abs(vec[5])+abs(vec[1]))/2)>=0.1){
            count <- count + 1
            if (count >= 3){
              vec <- vec[-5]
            }
          }
        }
      }
    }
  }
  
  if(length(vec)==6){
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[1])-abs(vec[5]))/((abs(vec[1])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[1])-abs(vec[6]))/((abs(vec[1])+abs(vec[6]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-1]
              }
            }
          }
        }
      }
    }
    if (count >= 4){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[2])-abs(vec[5]))/((abs(vec[2])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[2])-abs(vec[6]))/((abs(vec[2])+abs(vec[6]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-2]
              }
            }
          }
        }
      }
    }
    if (count >= 4){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[3])-abs(vec[5]))/((abs(vec[3])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[3])-abs(vec[6]))/((abs(vec[3])+abs(vec[6]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-3]
              }
            }
          }
        }
      }
    }
    if (count >= 4){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[4])-abs(vec[5]))/((abs(vec[4])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[4])-abs(vec[6]))/((abs(vec[4])+abs(vec[6]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-4]
              }
            }
          }
        }
      }
    }
    if (count >= 4){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[5])-abs(vec[2]))/((abs(vec[5])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[5])-abs(vec[3]))/((abs(vec[5])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[5])-abs(vec[4]))/((abs(vec[5])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[5])-abs(vec[1]))/((abs(vec[5])+abs(vec[1]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[5])-abs(vec[6]))/((abs(vec[5])+abs(vec[6]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-5]
              }
            }
          }
        }
      }
    }
    if (count >= 4){
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[6])-abs(vec[2]))/((abs(vec[6])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
      if((abs(vec[6])-abs(vec[3]))/((abs(vec[6])+abs(vec[3]))/2)>=0.1){
        count <- count + 1
        if((abs(vec[6])-abs(vec[4]))/((abs(vec[6])+abs(vec[4]))/2)>=0.1){
          count <- count + 1
          if((abs(vec[6])-abs(vec[5]))/((abs(vec[6])+abs(vec[5]))/2)>=0.1){
            count <- count + 1
            if((abs(vec[6])-abs(vec[1]))/((abs(vec[6])+abs(vec[1]))/2)>=0.1){
              count <- count + 1
              if (count >= 4){
                vec <- vec[-6]
              }
            }
          }
        }
      }
    }
  }
  
  if(length(vec)==7){
    if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[5]))/((abs(vec[1])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[6]))/((abs(vec[1])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[7]))/((abs(vec[1])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-1]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[5]))/((abs(vec[2])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[6]))/((abs(vec[2])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[7]))/((abs(vec[2])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-2]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[5]))/((abs(vec[3])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[6]))/((abs(vec[3])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[7]))/((abs(vec[3])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-3]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[5]))/((abs(vec[4])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[6]))/((abs(vec[4])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[7]))/((abs(vec[4])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-4]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[5])-abs(vec[2]))/((abs(vec[5])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[3]))/((abs(vec[5])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[4]))/((abs(vec[5])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[1]))/((abs(vec[5])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[6]))/((abs(vec[5])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[7]))/((abs(vec[5])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-5]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[6])-abs(vec[2]))/((abs(vec[6])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[3]))/((abs(vec[6])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[4]))/((abs(vec[6])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[5]))/((abs(vec[6])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[1]))/((abs(vec[6])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[7]))/((abs(vec[6])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-6]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[7])-abs(vec[2]))/((abs(vec[7])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[3]))/((abs(vec[7])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[4]))/((abs(vec[7])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[5]))/((abs(vec[7])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[6]))/((abs(vec[7])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[1]))/((abs(vec[7])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 4){
      vec <- vec[-7]
    }
  }
  
  if(length(vec)==8){
    if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[5]))/((abs(vec[1])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[6]))/((abs(vec[1])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[7]))/((abs(vec[1])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[8]))/((abs(vec[1])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-1]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[5]))/((abs(vec[2])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[6]))/((abs(vec[2])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[7]))/((abs(vec[2])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[8]))/((abs(vec[2])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-2]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[5]))/((abs(vec[3])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[6]))/((abs(vec[3])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[7]))/((abs(vec[3])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[8]))/((abs(vec[3])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-3]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[5]))/((abs(vec[4])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[6]))/((abs(vec[4])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[7]))/((abs(vec[4])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[8]))/((abs(vec[4])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-4]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[5])-abs(vec[2]))/((abs(vec[5])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[3]))/((abs(vec[5])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[4]))/((abs(vec[5])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[1]))/((abs(vec[5])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[6]))/((abs(vec[5])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[7]))/((abs(vec[5])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[8]))/((abs(vec[5])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-5]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[6])-abs(vec[2]))/((abs(vec[6])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[3]))/((abs(vec[6])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[4]))/((abs(vec[6])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[5]))/((abs(vec[6])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[1]))/((abs(vec[6])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[7]))/((abs(vec[6])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[8]))/((abs(vec[6])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-6]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[7])-abs(vec[2]))/((abs(vec[7])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[3]))/((abs(vec[7])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[4]))/((abs(vec[7])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[5]))/((abs(vec[7])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[6]))/((abs(vec[7])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[1]))/((abs(vec[7])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[8]))/((abs(vec[7])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-7]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[8])-abs(vec[2]))/((abs(vec[8])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[3]))/((abs(vec[8])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[4]))/((abs(vec[8])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[5]))/((abs(vec[8])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[6]))/((abs(vec[8])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[7]))/((abs(vec[8])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[1]))/((abs(vec[8])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-8]
    }
  }
      
  if(length(vec)==9){
    if((abs(vec[1])-abs(vec[9]))/((abs(vec[1])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[2]))/((abs(vec[1])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[3]))/((abs(vec[1])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[4]))/((abs(vec[1])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[5]))/((abs(vec[1])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[6]))/((abs(vec[1])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[7]))/((abs(vec[1])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[1])-abs(vec[8]))/((abs(vec[1])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-1]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[2])-abs(vec[1]))/((abs(vec[2])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[9]))/((abs(vec[2])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[3]))/((abs(vec[2])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[4]))/((abs(vec[2])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[5]))/((abs(vec[2])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[6]))/((abs(vec[2])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[7]))/((abs(vec[2])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[2])-abs(vec[8]))/((abs(vec[2])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-2]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[3])-abs(vec[1]))/((abs(vec[3])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[9]))/((abs(vec[3])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[2]))/((abs(vec[3])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[4]))/((abs(vec[3])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[5]))/((abs(vec[3])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[6]))/((abs(vec[3])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[7]))/((abs(vec[3])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[3])-abs(vec[8]))/((abs(vec[3])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-3]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[4])-abs(vec[1]))/((abs(vec[4])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[9]))/((abs(vec[4])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[2]))/((abs(vec[4])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[3]))/((abs(vec[4])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[5]))/((abs(vec[4])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[6]))/((abs(vec[4])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[7]))/((abs(vec[4])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[4])-abs(vec[8]))/((abs(vec[4])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-4]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[5])-abs(vec[1]))/((abs(vec[5])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[9]))/((abs(vec[5])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[2]))/((abs(vec[5])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[4]))/((abs(vec[5])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[3]))/((abs(vec[5])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[6]))/((abs(vec[5])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[7]))/((abs(vec[5])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[5])-abs(vec[8]))/((abs(vec[5])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-5]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[6])-abs(vec[1]))/((abs(vec[6])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[9]))/((abs(vec[6])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[2]))/((abs(vec[6])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[4]))/((abs(vec[6])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[5]))/((abs(vec[6])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[3]))/((abs(vec[6])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[7]))/((abs(vec[6])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[6])-abs(vec[8]))/((abs(vec[6])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-6]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[7])-abs(vec[1]))/((abs(vec[7])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[9]))/((abs(vec[7])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[2]))/((abs(vec[7])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[4]))/((abs(vec[7])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[5]))/((abs(vec[7])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[6]))/((abs(vec[7])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[3]))/((abs(vec[7])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[7])-abs(vec[8]))/((abs(vec[7])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-7]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[8])-abs(vec[1]))/((abs(vec[8])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[9]))/((abs(vec[8])+abs(vec[9]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[2]))/((abs(vec[8])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[4]))/((abs(vec[8])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[5]))/((abs(vec[8])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[6]))/((abs(vec[8])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[7]))/((abs(vec[8])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[8])-abs(vec[8]))/((abs(vec[8])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-8]
      return(vec)
    }
    else{
      count <- 0
    }
    if((abs(vec[9])-abs(vec[1]))/((abs(vec[9])+abs(vec[1]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[2]))/((abs(vec[9])+abs(vec[2]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[3]))/((abs(vec[9])+abs(vec[3]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[4]))/((abs(vec[9])+abs(vec[4]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[5]))/((abs(vec[9])+abs(vec[5]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[6]))/((abs(vec[9])+abs(vec[6]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[7]))/((abs(vec[9])+abs(vec[7]))/2)>=0.1){
      count <- count + 1
    }
    if((abs(vec[9])-abs(vec[8]))/((abs(vec[9])+abs(vec[8]))/2)>=0.1){
      count <- count + 1
    }
    if (count >= 5){
      vec <- vec[-9]
    }
  }
  vec
}

row.names(result) <- NULL

#Find gene of intrest index
m<-1
n<-2
while(identical(result$gene[m], result$gene[n])){
  m<-m+1
  n<-n+1
  if (m > nrow(result)||n > nrow(result)){
    break
  }
}

#Calculate means and omit outliers
x <- m + 1
test <- x
for (i in 1:nrow(result)){
  test <- test + 1
  if (test > nrow(result)){
    break
  }
  while(result$technical_replicate[x] != 1){
    x <- x + 1
    if (x > nrow(result)){
      break
    }
  }
  x <- x + 1
  if(result$technical_replicate[x] == 1){
    x <- x - 1
    result$mean[x] <- result$Dct[x]
    x <- x + 1
    y <- y + 1
    next
  }
  x <- x - 1
  if (result$technical_replicate[test] == 1){
    test = test+1
  }
  while(result$technical_replicate[test] != 1){
    y <- test
    test <- test + 1
    if (test > nrow(result)){
      y <- test - 1
      break
    }
  }
  if ((x > nrow(result))||(y > nrow(result))){
    break
  }
  if (result$Dct[x] == 0 || result$Dct[y] == 0){
    x <- x + 1
    y <- y + 1
    next
  }
  if (omit == 1){
    cor1 <- find.out(result$Dct[x:y])
    for (i in 1:length(cor1)) {
      cor1 <- find.out(cor1)
    }
    result$mean[y]<-mean(cor1)
  }
  else if (omit == 0){
    result$mean[y] <- mean(result$Dct[x:y])
  }
  x <- x+1
}

#Save Delta CT
direc <- choose.dir(default = "", caption = "Select folder to save data")
for(z in 1:nrow(gene)){
  for (y in 1:nrow(result)){
    if (result$gene[y] == gene$Gene.1[z]){
      marker1 <- y
      while (result$gene[y] == gene$Gene.1[z]) {
        y <- y + 1
        if(y > nrow(result)){
          break
        }
      }
      marker2 <- y - 1
    }
    if (exists("marker1") || exists("marker2")){
      break
    }
  }
  for (x in 1:nrow(result)){
    if (result$gene[x] == gene$Gene.2[z]){
      marker3 <- x
      while (result$gene[x] == gene$Gene.2[z]) {
        x <- x + 1
        if(x>nrow(result)){
          break
        }
      }
      marker4 <- x - 1
    }
    if (exists("marker3") || exists("marker4")){
      break
    }
  }
  write.table(result[c(marker1:marker2,marker3:marker4),], paste0(direc,  "\\", gene[z,2], "_DCT.txt"), sep="\t",row.names=FALSE)
  rm(marker1)
  rm(marker2)
  rm(marker3)
  rm(marker4)
  write.table(result, paste0(direc,  "\\DCT.txt"), sep="\t",row.names=FALSE)
}

#Display results
View(result)