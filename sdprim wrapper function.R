# to install old version of sdtoolkit: remotes::install_url(url="https://cran.r-project.org/src/contrib/Archive/sdtoolkit/sdtoolkit_2.31.tar.gz")

sdprim


library(sdtoolkit)
# library(ecr)
library(nsga2R)

prim=function(x,y, thresh=NULL, box.init = NULL, peel.alpha = .05, paste.alpha =.05, mass.min = 0.01,
              threshold = 0, pasting = TRUE, verbose = F,
              threshold.type = 1, paste.all = T, coverage = TRUE,
              showbounds = TRUE, style = "ineq", npts = NA, ninter = NA,
              nbump = 10, repro = F, dfrac = 0.5, peel_crit=2){
  
  
  ################ prim.traj() in non-function form #################
  
  options(digits = 4)
  d <- ncol(x)
  n <- nrow(x)
  k.max <- ceiling(1/mass.min)
  num.boxes <- k.max
  y.mean <- mean(y)
  mass.init <- length(y)/n
  
  if (is.null(box.init)) { 

      box.init <- apply(x, 2, range)
      box.diff <- box.init[2, ] - box.init[1, ]
      box.init[1, ] <- box.init[1, ] - 10 * paste.alpha * 
      box.diff
      box.init[2, ] <- box.init[2, ] + 10 * paste.alpha * 
      box.diff
      
    # }
    

  }
  
  ############## end prim.traj () in non-function form ##################
  
  # find.traj() returns the boxes and there bounds, with no other info
  boxseq <- find.traj(x = data.frame(x), y = y, box = box.init, peel.alpha = peel.alpha, 
                      paste.alpha = paste.alpha, mass.min = mass.min, threshold = threshold, d = d, n = n, pasting = pasting, 
                      verbose = verbose, paste.all = paste.all, peel_crit = peel_crit)
  
  # min(y)-0.1 * abs(min(y))
  
  # traj.info() returns the boxe constraints, plus x and y data for each box, plus T / F indicating what dimensions are constrained
  # also provides y.mean (I think same as density), and coverage (relcoverage)
  trajinf <- traj.info(x = data.frame(x), y = y, box.seq = boxseq, 
                       npts = npts, ninter = ninter)
  
  # all the boxes I want more data for
  checpts=1:length(boxseq$box)
  
  # computes pvalues for each dimension of each box. Usually uses checpts interactively selected in trajplot(), but I want all data for ALL boxes in boxseq
  #, so i replaced checpts with 1:length(boxseq$box)
  pvallist <- pvallister(checpts, x, y, trajinf)
  
  boxdata=trajinf
  boxdata[["pvals"]]=pvallist

  return(boxdata)
      
}


multiobjective_prim=function(x, y, qp_val=.1, peel_crit=2){
  
  # run PRIM using ALL metrics
  # significant=c()
  # for (i in 1:ncol(x)){
  #   x_temp=matrix(x[,i])
  #   boxData=prim(x_temp,y, peel_crit = peel_crit)
  #   pvals=boxData$pvals
  #   
  #   if (TRUE %in% (pvals[[1]] < qp_val)){
  #     significant=c(significant, colnames(x)[i])
  #   }
  #   
  #   
  # }

  boxData=prim(x,y, peel_crit = peel_crit)
  pvals=boxData$pvals
  dimlist=boxData$dimlist
  dim_names=names(boxData$dimlist[[1]]$either) # get metric names

  # find which metrics have qp values meeting the threshold

  significant=c()

  for (i in 1:length(pvals)){ # loop through pval matrices
    
    if (TRUE %in% (pvals[[i]] < qp_val)){ # if there is 1 or more significant box limit, find the metric and save it

      rc=which(pvals[[i]] < qp_val, arr.ind = T) # get row and column indices
      r=rc[,1] # get just the row. The row corresponds to the dimension

      used_dims=dimlist[[i]]$either[which(dimlist[[i]]$either==TRUE)] # returns all dimensions used in the box, both significant and non-significant
      dimensions=names(used_dims[r]) # get the name of only the significant dimension

      significant=c(significant, dimensions)

    }

    significant=unique(significant) # use unique() to remove duplicates

  }
  # significant=c("Driest10yrAVG", "LowQ.MaxDur.hist")
  
  
  # use combn() to create all combinations of 1 to D significant dimensions
  combos=list()
  
  for (i in 1:length(significant)){
    
    combos[[i]]=combn(x=significant, m=i)
    
  }
  
  # run PRIM on all combinations
  
  comboPRIM=list()
  c=1

  for (d in 1:length(combos)){ # loop through different number of dimensions (1D box, 5D box, etc)
    
    for (i in 1:ncol(combos[[d]])){ # loop through different combos of dimensions for each different box size
      
      # subset x data to contain only current combinations of dimensions
      
      dims_temp=combos[[d]][,i]
      x_temp=data.frame(x[, c(dims_temp)])
      colnames(x_temp)=dims_temp
      
      prim_temp=prim(x=x_temp, y=y,peel_crit = peel_crit)
      pvals=prim_temp$pvals
      
      remove_boxes=c()
      
      # remove boxes not meeting qp threshold on BOTH DIMENSIONs
      for (p in 1:length(pvals)){ # check pvals loop
        
        TF=pvals[[p]] < qp_val # Convert to TF. False means NOT significant
        number_insignificant_dimensions=sum((rowSums(TF)<1)) # where the row sum is less than 1, none of the box constrainsts on that dimension are significant. By summing, you know how many dimensions are insignificant
        
        if (number_insignificant_dimensions > 0){
          
          remove_boxes=c(remove_boxes, p)

        } # end if number_insignificant_dimensions > 0
        
      } # end check pvals loop
      
      if (length(prim_temp$x) == length(remove_boxes)){ # to avoid throwing errors by trying to remove all elements, 
        # simply don't save prim_temp in comboPRIM by using next, which goes to a new combination of dimensions without
        # saving prim_temp in comboPRIM
        next
      } else { # remove boxes from prim_temp
        
        if (length(remove_boxes)>0){
          
          for(element in names(prim_temp)){ # list for loop
            
            if (typeof(prim_temp[[element]])=="double"){
              
              prim_temp[[element]]=prim_temp[[element]][-remove_boxes]
              
            } else {
              prim_temp[[element]][remove_boxes]=NULL
            }
            
          } # end list for loop
          
        }
 
        comboPRIM[[c]]=prim_temp
        c=c+1
        
      } # end remove boxes
      
    } # end loop through combos of dimensions for each different box size loop
    
  } # end loop through different numbers of dimensions
  
  #################### non dominated filter #####################
  
  ### put all boxes into matrix form with coverage, density, and number of dimensions
  
  coverage=c()
  density=c()
  dimensions=c()
  
  for(i in 1:length(comboPRIM)){ # loop through prim iterations stored in comboPRIM
    for (b in 1:length(comboPRIM[[i]]$y.mean)){ # loop through boxes
      coverage=c(coverage, comboPRIM[[i]]$relcoverage[b])
      density=c(density, comboPRIM[[i]]$y.mean[b])
      dimensions=c(dimensions, nrow(comboPRIM[[i]][["pvals"]][[b]]))
    }
  }
  
  temp.metrics=data.frame(coverage=coverage*(-1), density=density*(-1), dimensions) # *-1 to minimize coverage and density
  
  ### perform non-domination operator
  
  nd.list=fastNonDominatedSorting(temp.metrics)
  front1=nd.list[[1]]
  
  temp.metrics=temp.metrics[front1,]
  
  box.metrics=data.frame(coverage=temp.metrics$coverage*(-1), density=temp.metrics$density*(-1), dimensions=temp.metrics$dimensions)
  
  ### extract non-dominated boxes from comboPRIM
  
  ## first, put all boxes into one concatanated object
  x_temp=c()
  y_temp=c()
  y.mean=c()
  box=c()
  mass=c()
  dimlist=c()
  relcoverage=c()
  pvals=c()
  for (i in 1:length(comboPRIM)){
    x_temp=c(x_temp, comboPRIM[[i]]$x)
    y_temp=c(y_temp, comboPRIM[[i]]$y)
    y.mean=c(y.mean, comboPRIM[[i]]$y.mean)
    box=c(box, comboPRIM[[i]]$box)
    mass=c(mass, comboPRIM[[i]]$mass)
    dimlist=c(dimlist, comboPRIM[[i]]$dimlist)
    relcoverage=c(relcoverage, comboPRIM[[i]]$relcoverage)
    pvals=c(pvals, comboPRIM[[i]]$pvals)
  }

  comboPRIM_concat=list(x=x_temp, y=y_temp, y.mean=y.mean,
                        box=box, dimlist=dimlist, relcoverage=relcoverage, pvals=pvals )
  
  ## extract
  
  return_PRIM=list()
  for (element in names(comboPRIM_concat)){
  return_PRIM[[element]]=comboPRIM_concat[[element]][front1]
  }

  return(list(box.metrics=box.metrics, box.info=return_PRIM))
  
    
}

################################# defineCOI function ######################

defineCOI=function(objectives_names=input$COIObj,  objectives_data=obj$filter, thresholds=thresholds, AndOrTF=input$AndOrTF){
  
  calcs=matrix(ncol=length(objectives_names)+2, nrow=nrow(objectives_data))
  
  for (cr in 1:length(objectives_names)){
    
      calcs[,cr]=ifelse(objectives_data[,which(colnames(objectives_data)==objectives_names[cr])] > thresholds[cr], 1,0) # a 1 means Case of Interest (COI)

  } # objectives thresholds loop
  
  if (length(objectives_names)==1){ # only one objective, don't need to do rowSums operator
    calcs[,length(objectives_names)+1]=calcs[,1]
  } else {
    calcs[,length(objectives_names)+1]=rowSums(calcs[,1:length(objectives_names)])
  }
  
  
  if (AndOrTF == TRUE){ # AND
    calcs[,ncol(calcs)]=ifelse(calcs[,(ncol(calcs)-1)] == length(objectives_names), 1, 0) # if sum is equal to number of objectives, then union condition is true (1)
  } else { # OR
    calcs[,ncol(calcs)]=ifelse(calcs[,(ncol(calcs)-1)] > 0, 1, 0) # if sum is greater than 0, then intersection condition is TRUE (1)
  }
  
  COI=calcs[,ncol(calcs)]

  return(COI)
    
}









