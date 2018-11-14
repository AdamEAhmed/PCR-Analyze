#Clear
rm(list = ls())

##SELECT DCT FILE##
result <- read.table(choose.files(caption = "SELECT DCT FILE"), header = TRUE, sep="\t")

##SELECT GENE FILE##
gene <- read.table(choose.files(caption = "SELECT GENE FILE"), header = TRUE, sep="\t")

##SELECT GROUPS FILE##
instr <- read.table(choose.files(caption = "SELECT GROUPS FILE"), header = TRUE, sep="\t")

###Set to 1 to omit outliers or set to 0 to include them. (default 1)###
omit <- 1

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

#Extract means  and convert to numeric
extracted <- subset(result, result$mean!="Omitted", select = 8)
names(extracted)<-"Dct"
extracted$Dct <- as.numeric(as.character(extracted$Dct))
mer <- merge(x = result, y = extracted, by = 0, sort = FALSE)

final <- subset(mer, select = c(4:7,10))

#Sort groups
temp1 <- data.frame(matrix(ncol = 1, nrow = 0))
for(i in 1:length(instr)){
  temp2 <- data.frame(instr[,i])
  names(temp1) <- names(temp2)
  temp1 <- rbind(temp1, temp2)
}
final <- final[order(match(final$treatment, temp1[,])),]
row.names(final)<-NULL
final <- final[order(match(final$gene, gene[,2])),]
row.names(final)<-NULL

#Label biological replicates
names(final)[4]<-"biological_replicate"
m <- 1
n <- 2
for(i in 1:nrow(final)){
  if(m > nrow(final) || n > nrow(final)){
    break
  }
  count <- 1
  while(identical(final$treatment[m], final$treatment[n])){
    final$biological_replicate[m] <- count
    count <- count + 1
    m <- m + 1
    n <- n + 1
  }
  if(!identical(final$treatment[m], final$treatment[n])){
    final$biological_replicate[m] <- count
    m <- m + 1
    n <- n + 1
  }
}

if("Ignored" %in% colnames(instr)){
  testing <- data.frame(matrix(ncol = 1, nrow = nrow(final)))
  for (z in 1:nrow(final)) {
    testing[z,1] <- instr[1,length(instr)] == final[z,3]
  }
  tempFinal <- final[testing[,1],]
  final <- final[!testing[,1],]
}

#Average of means
statistics<-setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Gene", "Treatment", "DDct", "2^DDCt", "SEM", "UCL", "LCL"))
if(omit == 1){
  final$omitted<-NA
}
final$avg_mean <- NA
i <- 1
j <- 2
a <- 1
for(k in 1:nrow(final)){
  intial <- i
  while(identical(final$treatment[i],final$treatment[j])) {
    i <- i + 1
    j <- j + 1
    if (i > nrow(final)|| j > nrow(final)){
      break
    }
  }
  if(omit == 1){
    cor2 <- find.out(final$Dct.y[intial:i])
    for (l in 1:length(cor2)) {
      cor2 <- find.out(cor2)
    }
    final$omitted[intial:i] <- !final$Dct.y[intial:i] %in% cor2
    final$avg_mean[i]<-mean(cor2)
    final$sd[i] <- sd(cor2)
    final$count[i] <- length(cor2)
  }
  else if(omit == 0){
    final$avg_mean[i] <- mean(final$Dct.y[intial:i])
    final$sd <- sd(final$Dct.y[intial:i])
    final$count[i] <- length(final$Dct.y[intial:i])
  }
  if(i>= nrow(final)){
    break
  }
  i <- i + 1
  j <- j + 1
}

#Calculate delta-delta ct
k <- 1
for (z in 1:nrow(gene)) {
  ind <- gene[z,3]
  for (a in 1:nrow(final)) {
    if (final$gene[a]==gene$Gene.2[z]){
      for (b in a:nrow(final)){
        if(final$treatment[b]==instr[1,gene[z,3]]){
          c <- b
          while (is.na(final$avg_mean[c])){
            c <- c + 1
            if(c > nrow(final)){
              break
            }
          }
          control <- final$avg_mean[c]
          sdc <- final$sd[c]
          i <- 1
          for (l in 1:nrow(instr)) {
            for (g in 1:nrow(final)) {
              if(i>nrow(final)){
                break
              }
              if (final$treatment[i]==instr[l,ind]){
                while (is.na(final$avg_mean[i])){
                  i <- i + 1          
                }
                statistics[k,1] <- paste0(final$gene[i])
                statistics[k,2] <- paste0(final$treatment[i])
                statistics[k,3] <- final$avg_mean[i] - control
                statistics[k,5] <- sqrt(final$sd[i]^2+sdc^2)/sqrt(final$count[i])
                i<-i+1
                k <- k + 1
                break
              }
              else{
                i<-i+1
              }
            }
          }
          break
        }
      }
    }
    break
  }
}

for (i in 1:nrow(statistics)) {
  if (statistics$DDct[i] == 0){
    statistics$SEM[i] <- 0
  }
}

#Calculate 2^DDCt, UCL, & LCL
for (i in 1:nrow(statistics)) {
  statistics[i,4] <- 2^(statistics$DDct[i])
  statistics$LCL[i] <- 2^(statistics$DDct[i]-(1.96*statistics$SEM[i]))
  statistics$UCL[i] <- 2^(statistics$DDct[i]+(1.96*statistics$SEM[i]))
  if(statistics$Treatment[i] %in% instr[1,] && statistics$DDct[i]==0){
    statistics$SEM[i] <- 0
  }
}

fix.labels <- function(x){
  x+1
}

direc <- choose.dir(default = "", caption = "Select folder to save plots and data")

instr <- cbind(instr, "index"=1:nrow(instr))

#Plot Graph
library(ggplot2)
m <- 1
n <- 2
for (i in 1:nrow(gene)) {
  intial <- m
  while (statistics$DDct[n]!=0){
    n <- n + 1
  }
  m <- n
  n <- n - 1
  if(min(statistics[intial:n,3:7])>-5 && max(statistics[intial:n,3:7])<5){
    scale <- 0.1
  }
  else if(min(statistics[intial:n,3:7])>-14 && max(statistics[intial:n,3:7])<14){
    scale <- 0.5
  }
  else{
    scale <- 1.0
  }
  Plot <- ggplot(statistics[intial:n,3:7]) +
    xlab("Group") +
    ylab(paste0("-Fold Change ", gene[i,1], "/",  gene[i,2]))+
    geom_bar( aes(x=reorder(instr[1:length(intial:n),gene[i,3]], instr$index[1:length(intial:n)]), y=statistics$"2^DDCt"[intial:n]-1), stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=1:nrow(statistics[intial:n,3:7]), ymin=statistics$LCL[intial:n]-1, ymax=statistics$UCL[intial:n]-1), width=0.4, colour="orange", alpha=0.9, size=1.5) +
    scale_y_continuous(breaks = (seq((round(min(statistics[intial:n,3:7]), digits = 1)-1), (round(max(statistics[intial:n,3:7]), digits = 1)+1), by = scale)), labels=fix.labels)
  ##SAVE GRAPH##
  ggsave(paste0(direc,  "\\", gene[i,1], "_", gene[i,2], "-", instr[1,gene[i,3]], "_plot.jpg"), plot=Plot)
  n <- n + 1
}

if(exists("tempFinal")){
  final <- merge(final, tempFinal, all = TRUE)
  final <- final[order(match(final$treatment, temp1[,])),]
  row.names(final)<-NULL
  final <- final[order(match(final$gene, gene[,2])),]
  row.names(final)<-NULL
}

##Save Data
m <- 1
n <- 2
for (i in 1:nrow(gene)) {
  intial <- m
  while (statistics$DDct[n]!=0){
    n <- n + 1
  }
  m <- n
  n <- n - 1
  write.table(statistics[intial:n,], paste0(direc,  "\\", gene[i,2], "-", instr[1,gene[i,3]], "_stats.txt"), sep="\t",row.names=FALSE)
  n <- n + 1
}
write.table(final, paste0(direc,  "\\", "Final.txt"), sep="\t", row.names=FALSE)

#Display results
View(final)
View(statistics)