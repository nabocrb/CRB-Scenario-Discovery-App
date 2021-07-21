library(nsga2R)
library(dplyr)
library(ggplot2)
library(reshape2)
library(plotly)
library(shiny)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(GGally)
library(grid)
library(ggpubr)
library(cowplot)
library(rlang) # as_string
library(stringr)
library(forcats) # to order violin plots according to user input
library(ggnewscale)

if(exists("trigger")){rm("trigger")}

############################## data ###############################

### load stress testing data from robustness runs

obj_all=read.table('objectives_all463.txt') # load data for calculating robustness
obj_all=rename(obj_all, "ID"="policy")

vol_index=which(colnames(obj_all) %in% c('LB.Shortage.Volume',"LB.Shortage.Volume.Policy", "Max.Annual.LB.Shortage"))
obj_all[,vol_index]=obj_all[,vol_index]/1000 # convert to KAF
obj_all[["Powell.WY.Release"]]=obj_all[["Powell.WY.Release"]]/(10^6) # convert to MAF

### load uncertainty metrics

SOW_500_300_100=readRDS("SOW ensemble 500_300_100.rds") # list of SOW ensembles, for sensitivity analysis of 300 vs 500 ensemble
SOW_id=sort(SOW_500_300_100[[2]]$model) # SOW IDs of the 300 member SOW ensemble

############################## pquantile #########################

pquantile=function(x, alpha, ptype){
  
  N <- length(x)
  xs <- sort(x)
  intm <- alpha * N
  pthresh <- 0.5 * (xs[max(floor(intm), 1)] + xs[max(ceiling(intm), 
                                                     1)])
  if (xs[max(floor(intm), 1)] == xs[max(ceiling(intm), 1)] & 
      min(xs) != max(xs)) {
    if ((floor(intm) <= 1 | ceiling(intm) >= N)) {
      if (ptype == "lowend") {
        i <- 2
        while ((i <= N) & (xs[i] == xs[i - 1])) {
          i <- i + 1
        }
        if (i <= N) {
          pthresh <- 0.5 * (xs[i] + xs[i - 1])
        }
        else {
          pthresh <- NA
        }
      }
      if (ptype == "highend") {
        i <- N - 1
        while ((i > 1) & (xs[i] == xs[i + 1])) {
          i <- i - 1
        }
        if (i > 0) {
          pthresh <- 0.5 * (xs[i] + xs[i + 1])
        }
        else {
          pthresh <- NA
        }
      }
    }
    else if (ptype == "lowend") {
      if (xs[floor(intm)] == xs[floor(intm) - 1]) {
        i <- ceiling(intm)
        while ((i <= N) & (xs[i] == xs[i - 1])) {
          i <- i + 1
        }
        if (i <= N) {
          pthresh <- 0.5 * (xs[i] + xs[i - 1])
        }
        else {
          pthresh <- NA
        }
      }
    }
    else if (ptype == "highend") {
      if (xs[floor(intm)] == xs[floor(intm) + 1]) {
        i <- floor(intm)
        while ((i > 1) & (xs[i] == xs[i + 1])) {
          i <- i - 1
        }
        if (i > 0) {
          pthresh <- 0.5 * (xs[i] + xs[i + 1])
        }
        else {
          pthresh <- NA
        }
      }
    }
  }
  return(pthresh)

}

##################### vol.box ##############################

vol.box=function (box, x){
  x.ind.curr <- rep(TRUE, nrow(x))
  box.curr <- box
  for (j in 1:ncol(x)) {
    x.ind.curr <- x.ind.curr & (x[, j] >= box.curr[1, j]) & 
      (x[, j] <= box.curr[2, j])
  }
  return(sum(x.ind.curr))
}


###################### peel.one #############################

peel.one=function (x, y, box, peel.alpha, mass.min, threshold, d, n, 
                   peel_crit){
  box.new <- box
  mass <- length(y)/n
  peel.alpha <- max(peel.alpha, 1/nrow(x))
  if (is.vector(x)){x=data.frame(x)}
    #return(NULL)
  y.mean <- mean(y)
  y.mean.peel <- matrix(0, nrow = 2, ncol = d)
  y.mean.peeled <- matrix(0, nrow = 2, ncol = d)
  box.vol.peel <- matrix(NA, nrow = 2, ncol = d)
  box.vol.peeled <- matrix(0, nrow = 2, ncol = d)
  box.supp.peeled <- matrix(0, nrow = 2, ncol = d)
  ranges <- apply(x, 2, range)
  spans <- ranges[2, ] - ranges[1, ]
  for (j in 1:d) {
    if (TRUE) {
      box.min.new <- pquantile(x[, j], peel.alpha, ptype = "lowend")
      box.max.new <- pquantile(x[, j], 1 - peel.alpha, 
                               ptype = "highend")
    }
    else {
      box.min.new <- min(x[, j])
      box.max.new <- max(x[, j])
      if (box.min.new != box.max.new) 
        stop("what?")
    }
    y.mean.peel[1, j] <- mean(y[x[, j] >= box.min.new])
    y.mean.peel[2, j] <- mean(y[x[, j] <= box.max.new])
    y.mean.peeled[1, j] <- mean(y[x[, j] < box.min.new])
    y.mean.peeled[2, j] <- mean(y[x[, j] > box.max.new])
    box.supp.peeled[1, j] <- sum(x[, j] < box.min.new)
    box.supp.peeled[2, j] <- sum(x[, j] > box.min.new)
    box.temp1 <- box
    box.temp2 <- box
    box.temp1.peeled <- box
    box.temp2.peeled <- box
    if (!is.na(box.min.new)) {
      box.temp1[1, j] <- box.min.new
      box.temp1.peeled[2, j] <- box.min.new
    }
    if (!is.na(box.max.new)) {
      box.temp2[2, j] <- box.max.new
      box.temp2.peeled[1, j] <- box.max.new
    }
    box.vol.peel[1, j] <- vol.box(box.temp1, x)
    box.vol.peel[2, j] <- vol.box(box.temp2, x)
    box.vol.peeled[1, j] <- vol.box(box.temp1.peeled, x)
    box.vol.peeled[2, j] <- vol.box(box.temp2.peeled, x)
  }
  toprint <- y.mean.peel
  if (peel_crit == 1) {
    y.mean.peel.max.ind <- which(y.mean.peel == max(y.mean.peel, 
                                                    na.rm = TRUE), arr.ind = TRUE)
  }
  else if (peel_crit == 2) {
    evaled.peel.crit2 <- (y.mean.peel - y.mean)/box.supp.peeled
    y.mean.peel.max.ind <- which(evaled.peel.crit2 == max(evaled.peel.crit2, 
                                                          na.rm = TRUE), arr.ind = TRUE)
  }
  else if (peel_crit == 3) {
    evaled.peel.crit3 <- (y.mean.peel - y.mean) * (nrow(x) - 
                                                     box.supp.peeled)/box.supp.peeled
    y.mean.peel.max.ind <- which(evaled.peel.crit3 == max(evaled.peel.crit3, 
                                                          na.rm = TRUE), arr.ind = TRUE)
  }
  else {
    stop("There is no implemented peeling criteria associated with the value \n\t\t\t\t\tthat was passed to the peel_crit argument")
  }
  nrr <- nrow(y.mean.peel.max.ind)
  if (nrr > 1) {
    box.vol.peel2 <- rep(0, nrr)
    for (j in 1:nrr) box.vol.peel2[j] <- box.vol.peel[y.mean.peel.max.ind[j, 
                                                                          1], y.mean.peel.max.ind[j, 2]]
    row.ind <- which(max(box.vol.peel2) == box.vol.peel2)[1]
  }
  else {
    row.ind <- 1
  }
  y.mean.peel.max.ind <- y.mean.peel.max.ind[row.ind, ]
  j.max <- y.mean.peel.max.ind[2]
  if (y.mean.peel.max.ind[1] == 1) {
    box.new[1, j.max] <- pquantile(x[, j.max], peel.alpha, 
                                   ptype = "lowend")
    x.index <- x[, j.max] >= box.new[1, j.max]
  }
  else if (y.mean.peel.max.ind[1] == 2) {
    box.new[2, j.max] <- pquantile(x[, j.max], 1 - peel.alpha, 
                                   ptype = "highend")
    x.index <- x[, j.max] <= box.new[2, j.max]
  }
  x.new <- x[x.index, ]
  y.new <- y[x.index]
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)
  if ((y.mean.new >= threshold) & (mass.new >= mass.min) & 
      (mass.new < mass) & (y.mean < 1)) 
    return(list(x = x.new, y = y.new, y.mean = y.mean.new, 
                box = box.new, mass = mass.new))
}

######################### in.box #################################

in.box=function (x, box, d, boolean = FALSE){
  x.box.ind <- rep(TRUE, nrow(x))
  for (i in 1:d) x.box.ind <- x.box.ind & (box[1, i] <= x[, 
                                                          i]) & (x[, i] <= box[2, i])
  if (boolean) 
    return(x.box.ind)
  else return(x[x.box.ind, ])
}


####################### paste.one ################################

paste.one=function (x, y, x.init, y.init, box, paste.alpha, mass.min, 
                    threshold, d, n){
  box.new <- box
  mass <- length(y)/n
  y.mean <- mean(y)
  n.box <- length(y)
  box.init <- apply(x.init, 2, range)
  if (is.vector(x)) 
    x <- as.matrix(t(x))
  y.mean.paste <- matrix(0, nrow = 2, ncol = d)
  mass.paste <- matrix(0, nrow = 2, ncol = d)
  box.paste <- matrix(0, nrow = 2, ncol = d)
  x.paste1.list <- list()
  x.paste2.list <- list()
  y.paste1.list <- list()
  y.paste2.list <- list()
  box.paste1 <- box
  box.paste2 <- box
  for (j in 1:d) {
    box.diff <- (box.init[2, ] - box.init[1, ])[j]
    box.paste1[1, j] <- box[1, j] - box.diff * paste.alpha
    box.paste2[2, j] <- box[2, j] + box.diff * paste.alpha
    x.paste1.ind <- in.box(x = x.init, box = box.paste1, 
                           d = d, boolean = TRUE)
    x.paste1 <- x.init[x.paste1.ind, ]
    y.paste1 <- y.init[x.paste1.ind]
    x.paste2.ind <- in.box(x = x.init, box = box.paste2, 
                           d = d, boolean = TRUE)
    x.paste2 <- x.init[x.paste2.ind, ]
    y.paste2 <- y.init[x.paste2.ind]
    if (box.diff > 0) {
      while (length(y.paste1) <= length(y) & box.paste1[1, 
                                                        j] >= box.init[1, j]) {
        box.paste1[1, j] <- box.paste1[1, j] - box.diff * 
          paste.alpha
        x.paste1.ind <- in.box(x = x.init, box = box.paste1, 
                               d = d, boolean = TRUE)
        x.paste1 <- x.init[x.paste1.ind, ]
        y.paste1 <- y.init[x.paste1.ind]
      }
      while (length(y.paste2) <= length(y) & box.paste2[2, 
                                                        j] <= box.init[2, j]) {
        box.paste2[2, j] <- box.paste2[2, j] + box.diff * 
          paste.alpha
        x.paste2.ind <- in.box(x = x.init, box = box.paste2, 
                               d = d, boolean = TRUE)
        x.paste2 <- x.init[x.paste2.ind, ]
        y.paste2 <- y.init[x.paste2.ind]
      }
    }
    y.mean.paste[1, j] <- mean(y.paste1)
    y.mean.paste[2, j] <- mean(y.paste2)
    mass.paste[1, j] <- length(y.paste1)/n
    mass.paste[2, j] <- length(y.paste2)/n
    x.paste1.list[[j]] <- x.paste1
    y.paste1.list[[j]] <- y.paste1
    x.paste2.list[[j]] <- x.paste2
    y.paste2.list[[j]] <- y.paste2
    box.paste[1, j] <- box.paste1[1, j]
    box.paste[2, j] <- box.paste2[2, j]
  }
  y.mean.paste.max <- which(y.mean.paste == max(y.mean.paste, 
                                                na.rm = TRUE), arr.ind = TRUE)
  if (nrow(y.mean.paste.max) > 1) {
    y.mean.paste.max <- cbind(y.mean.paste.max, mass.paste[y.mean.paste.max])
    y.mean.paste.max.ind <- y.mean.paste.max[order(y.mean.paste.max[, 
                                                                    3], decreasing = TRUE), ][1, 1:2]
  }
  else {
    y.mean.paste.max.ind <- as.vector(y.mean.paste.max)
  }
  j.max <- y.mean.paste.max.ind[2]
  if (y.mean.paste.max.ind[1] == 1) {
    x.new <- x.paste1.list[[j.max]]
    y.new <- y.paste1.list[[j.max]]
    box.new[1, j.max] <- box.paste[1, j.max]
  }
  else if (y.mean.paste.max.ind[1] == 2) {
    x.new <- x.paste2.list[[j.max]]
    y.new <- y.paste2.list[[j.max]]
    box.new[2, j.max] <- box.paste[2, j.max]
  }
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)
  if ((y.mean.new > threshold) & (mass.new >= mass.min) & 
      (y.mean.new >= y.mean) & (mass.new > mass)) 
    return(list(x = x.new, y = y.new, y.mean = y.mean.new, 
                box = box.new, mass = mass.new))
}

######################### find.traj ###############################


find.traj=function (x, y, box, peel.alpha, paste.alpha, mass.min, threshold, 
            d, n, pasting, verbose = FALSE, paste.all = FALSE, peel_crit){
    peel.traj <- list()
    y.mean <- mean(y)
    mass <- length(y)/n
    if ((y.mean >= threshold) & (mass >= mass.min)) 
      boxk.peel <- peel.one(x = x, y = y, box = box, peel.alpha = peel.alpha, 
                            mass.min = mass.min, threshold = threshold, d = d, 
                            n = n, peel_crit = peel_crit)
    else boxk.peel <- NULL
    boxk.temp <- NULL
    bi <- 0
    while (!is.null(boxk.peel)) {
      bi <- bi + 1
      boxk.temp <- boxk.peel
      peel.traj[[bi]] <- boxk.temp
      flush.console()
      boxk.peel <- peel.one(x = boxk.temp$x, y = boxk.temp$y, 
                            box = boxk.temp$box, peel.alpha = peel.alpha, mass.min = mass.min, 
                            threshold = threshold, d = d, n = n, peel_crit = peel_crit)
    }
    if (verbose) {
      cat("Peeling completed \n")
    }
    if (!pasting) {
      paste.traj <- lapply(peel.traj, "[[", "box")
    }
    else {
      paste.traj <- list()
      for (p in 1:length(peel.traj)) {
        boxk.paste <- peel.traj[[p]]
        while (!is.null(boxk.paste)) {
          boxk.temp <- boxk.paste
          boxk.paste <- paste.one(x = boxk.temp$x, y = boxk.temp$y, 
                                  box = boxk.temp$box, x.init = x, y.init = y, 
                                  paste.alpha = paste.alpha, mass.min = mass.min, 
                                  threshold = threshold, d = d, n = n)
        }
        if (verbose) {
          cat("Pasting completed\n")
        }
        paste.traj[[p]] <- boxk.temp$box
      }
    }
    boxk <- boxk.temp
    paste.seq <- list()
    paste.seq$box <- paste.traj
    paste.seq$num.class <- length(paste.seq$box)
    return(paste.seq)
}
  
###################### dimchecker ########################

dimchecker=function (x, box){
  box.init <- apply(x, 2, range)
  lows <- box[1, ] > box.init[1, ]
  highs <- box[2, ] < box.init[2, ]
  rdims <- (lows | highs)
  return(list(either = rdims, lower = lows, upper = highs))
}


####################### traj.info #########################
  
traj.info=function (x, y, box.seq, npts = NA, ninter = NA){
  m <- box.seq$num.class
  d <- ncol(x)
  n <- nrow(x)
  x.ind <- rep(TRUE, n)
  xy.list <- list()
  for (k in 1:m) {
    x.ind.curr <- x.ind
    box.curr <- box.seq$box[[k]]
    for (j in 1:d) {
      x.ind.curr <- x.ind.curr & (x[, j] >= box.curr[1, 
                                                     j]) & (x[, j] <= box.curr[2, j])
    }
    x.curr <- x[x.ind.curr & x.ind, ]
    box.mass.curr <- sum(x.ind.curr)/n
    xy.list$x[[k]] <- x.curr
    if (!missing(y)) {
      y.curr <- y[x.ind.curr & x.ind]
      y.mean.curr <- mean(y.curr)
      xy.list$y[[k]] <- y.curr
      xy.list$y.mean[[k]] <- y.mean.curr
    }
    xy.list$box[[k]] <- box.curr
    xy.list$mass[[k]] <- box.mass.curr
    xy.list$dimlist[[k]] <- dimchecker(x, box.curr)
  }
  precov <- xy.list$mass * xy.list$y.mean
  xy.list$relcoverage <- precov * nrow(x)/sum(y)
  xy.list$marcoverage <- precov * nrow(x)/ninter
  return(xy.list)
}
  
######################### boxconverter ####################

boxconverter=function (sdbox){
  dimvect <- c(1:length(sdbox$dimlist$either))[sdbox$dimlist$either]
  dimmat <- matrix(Inf, nrow = sum(sdbox$dimlist$either), 
                   ncol = 2)
  dimmat[, 1] <- -Inf
  
  test_var=length(sdbox$dimlist$lower)

  trigger=0
  
        for (i in 1:length(dimvect)) {
          
          test_var=sdbox$dimlist$lower[dimvect[i]]
          #cat(file=stderr(), 'test_var:', test_var, "\n") # prints dimension to R console
          
          if(!(test_var %in% c(TRUE, FALSE))){
            
            #showNotification("No boxes found. Try different PRIM settings or observable conditions.", type = "error", duration = NULL )
            
            trigger=1
            #cat(file=stderr(), 'function:', "boxconverter", "\n") # prints dimension to R console
            break
            
            }
            
            if (sdbox$dimlist$lower[dimvect[i]]) {
              dimmat[i, 1] <- sdbox$box[1, dimvect[i]]
            }
            if (sdbox$dimlist$upper[dimvect[i]]) {
              dimmat[i, 2] <- sdbox$box[2, dimvect[i]]
            }
    

        }
  
    if(trigger==1){
      return(NULL)
    } else {
      newbox <- list(dimvect, dimmat)
      return(newbox)
    }

  
}

######################### nullprob #########################

nullprob=function (dset, y = NULL, lbox){
  if (is.null(y)) {
    y <- dset[, ncol(dset)]
  }
  if (length(lbox[[1]]) > 1) {
    lvouts <- lvout(dset, y, lbox)
  }
  else {
    lvouts <- matrix(c(NA, NA, sum(y)/nrow(dset)), ncol = 3)
    vect <- dset[, lbox[[1]][1]]
    invect <- (vect > lbox[[2]][1, 1]) & (vect < lbox[[2]][1, 
                                                           2])
    attr(lvouts, "origtotin") <- sum(invect)
    attr(lvouts, "orighighin") <- sum(y[invect])
  }
  pvals <- vector(length = nrow(lvouts))
  for (i in 1:nrow(lvouts)) {
    pbase <- lvouts[i, 3]
    nbig <- attr(lvouts, "origtotin")
    intot <- attr(lvouts, "orighighin")
    ptot <- 0
    pvals[i] <- pbinom(intot - 1, nbig, pbase, lower.tail = FALSE)
  }
  return(pvals)
}

######################### lvout ############################

lvout=function (dset, y = NULL, lbox){
  rmdsizes <- vector(length = length(lbox[[1]]))
  rmdhighs <- vector(length = length(lbox[[1]]))
  for (i in 1:length(lbox[[1]])) {
    tempbox <- list()
    tempbox[[1]] <- lbox[[1]][-i]
    tempbox[[2]] <- lbox[[2]][-i, ]
    if (is.null(y)) {
      allin <- unionpts(dset, list(tempbox))
    }
    else {
      allin <- unionpts(cbind(dset, y), list(tempbox))
    }
    rmdsizes[i] <- sum(allin)
    if (is.null(y)) {
      rmdhighs[i] <- sum(dset[allin, ncol(dset)])
    }
    else {
      rmdhighs[i] <- sum(y[allin])
    }
  }
  rstats <- cbind(rmdsizes, rmdhighs, rmdhighs/rmdsizes)
  tempbox <- lbox
  if (is.null(y)) {
    allin <- unionpts(dset, list(tempbox))
  }
  else {
    allin <- unionpts(cbind(dset, y), list(tempbox))
  }
  origtotin <- sum(allin)
  if (is.null(y)) {
    orighighin <- sum(dset[allin, ncol(dset)])
  }
  else {
    orighighin <- sum(y[allin])
  }
  attr(rstats, "origtotin") <- origtotin
  attr(rstats, "orighighin") <- orighighin
  attr(rstats, "origdens") <- orighighin/origtotin
  return(rstats)
}

########################## unionpts ########################

unionpts=function (dset, blist){
  outdim <- ncol(dset)
  tpts <- nrow(dset)
  this <- sum(dset[, outdim])
  B <- length(blist)
  bdims <- rep(0, B)
  for (i in 1:B) {
    bdims[i] <- length(blist[[i]][[1]])
  }
  bincvecs <- matrix(TRUE, nrow = tpts, ncol = B)
  dnet <- c()
  for (b in 1:B) {
    dimvect <- blist[[b]][[1]]
    dnet <- c(dnet, dimvect)
    bmat <- blist[[b]][[2]]
    if (is.vector(bmat)) {
      bmat <- matrix(bmat, ncol = 2)
    }
    incvecs <- !logical(length = tpts)
    for (i in 1:length(dimvect)) {
      di <- dimvect[i]
      incvecs <- incvecs & (dset[, di] >= bmat[i, 1])
      incvecs <- incvecs & (dset[, di] < bmat[i, 2])
    }
    bincvecs[, b] <- incvecs
  }
  masvecs <- logical(length = nrow(dset))
  for (b in 1:B) {
    masvecs <- masvecs | bincvecs[, b]
  }
  return(masvecs)
}


######################### pvallister #######################
  
pvallister=function (checpts, x, y, trajinf){
  nboxes <- length(checpts)
  pvallist <- list()
  ofboxlist <- list()
  for (boxind in checpts) {
    boxy <- list(box = trajinf$box[[boxind]], dimlist = trajinf$dimlist[[boxind]])
    
    error_check=boxconverter(boxy)
 
    if(is.null(error_check)){
      
      trigger=1
      
      break
    } else {
      ofboxlist[[boxind]] <- boxconverter(boxy)
    }
    
  }
  
  if(exists("trigger")){
    return(NULL)
    cat(file=stderr(), 'function:', "pvallister", "\n") # prints dimension to R console
    
  } else {
    
    for (i in checpts) {
      pvallist[[i]] <- cbind(ofboxlist[[i]][[1]], nullprob(cbind(x, 
                                                                 y), y = NULL, ofboxlist[[i]]))
    }
    return(pvallist)
   
  } 

}

######################### prim #############################
  
prim=function(x,y, thresh=NULL, box.init = NULL, peel.alpha = .05, paste.alpha =.05, mass.min = 0.1,
                threshold = 0, pasting = TRUE, verbose = F,
                threshold.type = 1, paste.all = T, coverage = TRUE,
                showbounds = TRUE, style = "ineq", npts = NA, ninter = NA,
                nbump = 10, repro = F, dfrac = 0.5, peel_crit=2, remove_dims=F){
    
    
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
      box.init[1, ] <- box.init[1, ] - 1 * paste.alpha * 
        box.diff
      box.init[2, ] <- box.init[2, ] + 1 * paste.alpha * 
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
    
    if(is.null(pvallist)){
      return(NULL)
      cat(file=stderr(), 'function:', "prim", "\n") # prints dimension to R console
      
    } else {
      
    
    
    boxdata=trajinf
    boxdata[["pvals"]]=pvallist
    
    if (remove_dims==T){ ## remove unused dimensions from x and box elements
      
      for (i in 1:length(boxdata$x)){
        
        keep_dims=boxdata$dimlist[[i]]$either
        
        if (is.null(ncol(boxdata$x[[i]]))){
          
        } else {
          boxdata$x[[i]]=boxdata$x[[i]][,keep_dims]
          boxdata$box[[i]]=boxdata$box[[i]][,keep_dims]
        }
        
        
      }
      
    } 

    
    return(boxdata)
    
    }
    
  }

########################## MO.PRIM ######################################
  
  MO.PRIM=function(x, y, qp_val=.1, peel_crit=2, peel_alpha=input$PeelAlpha){
    
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
    
    boxData=prim(x,y, peel_crit = peel_crit, peel.alpha = peel_alpha)
    
    if(is.null(boxData)){
      showNotification("An error has occurred. Please try different PRIM settings or observable conditions.", type="error")
      cat(file=stderr(), 'function:', "MO.PRIM", "\n") # prints dimension to R console
      
      return(NULL)
      
    } else {
    
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
        
        prim_temp=prim(x=x_temp, y=y,peel_crit = peel_crit, peel.alpha = peel_alpha)
        pvals=prim_temp$pvals
        
        remove_boxes=c()
        
        # remove boxes not meeting qp threshold on BOTH DIMENSIONs
        for (p in 1:length(pvals)){ # check pvals loop
          
          TF=pvals[[p]] < qp_val # Convert to TF. False means NOT significant
          if(is.null(dim(TF))){
            TF=matrix(TF)
          }
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
  

    for (i in 1:length(comboPRIM)){

      
      if(i==1){
        x_temp=comboPRIM[[i]]$x
        y_temp=comboPRIM[[i]]$y
        y.mean=comboPRIM[[i]]$y.mean
        box=comboPRIM[[i]]$box
        mass=comboPRIM[[i]]$mass
        dimlist=comboPRIM[[i]]$dimlist
        relcoverage= comboPRIM[[i]]$relcoverage
        pvals=comboPRIM[[i]]$pvals
      } else {
        x_temp=c(x_temp, comboPRIM[[i]]$x)
        y_temp=c(y_temp, comboPRIM[[i]]$y)
        y.mean=c(y.mean, comboPRIM[[i]]$y.mean)
        box=c(box, comboPRIM[[i]]$box)
        mass=c(mass, comboPRIM[[i]]$mass)
        dimlist=c(dimlist, comboPRIM[[i]]$dimlist)
        relcoverage=c(relcoverage, comboPRIM[[i]]$relcoverage)
        pvals=c(pvals, comboPRIM[[i]]$pvals)
      }
      
      
    }
    
    comboPRIM_concat=list(x=x_temp, y=y_temp, y.mean=y.mean,
                          box=box, dimlist=dimlist, relcoverage=relcoverage, pvals=pvals )
    
    ## extract
    
    temp.box.info=list()
    for (element in names(comboPRIM_concat)){
      temp.box.info[[element]]=comboPRIM_concat[[element]][front1]
    }
    
    ## remove unused dimensions from x and box elements
    for (i in 1:length(temp.box.info$x)){
      
      keep_dims=temp.box.info$dimlist[[i]]$either
      
      if (is.null(ncol(temp.box.info$x[[i]]))){
        
      } else {
        temp.box.info$x[[i]]=temp.box.info$x[[i]][,keep_dims]
        temp.box.info$box[[i]]=temp.box.info$box[[i]][,keep_dims]
      }
      
      
    }
    
    
    
    ## rearrange so box indices are logical. Arrange by increasing dimension and density
    
    
    box.metrics$oldRow=1:nrow(box.metrics) # assign an old ID, which you will use to rearrange temp.box.info
    box.metrics=dplyr::arrange(box.metrics, dimensions, density) # rearrange so row index makes sense. Rearrange by increasing dimension, and increasing density within equal dimensions
    oldRow=box.metrics$oldRow
    
    box.info=list()
    for (element in names(temp.box.info)){ # loop through list elements
      for (i in 1:length(temp.box.info$x)){ # loop through boxes
        box.info[[element]][[i]]=temp.box.info[[element]][[oldRow[i]]]
      }
      
    }
    
    showNotification("Plotting")
    
    return(list(box.metrics=box.metrics[,-4], box.info=box.info))
    
    } # end if else
    
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
  
  
####################################### plotting #################################
  ################################################################################
  ################################################################################
  
### DV_plot
  
  metrics_4app_n500=readRDS('tradeoff_dataframes.rds')
  optimization=metrics_4app_n500$optimization
  
  bar_plot_data=readRDS(file='data for stacked bar plot.rds')
  long_data=bar_plot_data$long_data
  wide_data=bar_plot_data$wide_data
  wide_data$ID=1:nrow(wide_data)  

  ############################# stacked histogram Decison Variable plotting ###############################
  ##########################################################################################################
  
  # function to change ggplot legend size: https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
  addSmallLegend <- function(myPlot=bar_plot, pointSize = 0.25, textSize = 8, spaceLegend = 0.03) {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_blank(), #element_text(size = textSize)
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }
  
  DV_plot=function(long.data=long_data, wide.data=wide_data, metric= optimization,
                   to_plot=input$policyID, metric_label='order', preferred_direction='min', y_axis2=F,
                   labelsize=3, show_legend=F, height=275){ #labelsize=delayFontSlider()
    
    # filter for chosen policies
    # filter.long=dplyr::filter(long.data, policy %in% to_plot)
    # filter.wide=dplyr::filter(wide.data, ID %in% to_plot)
    # filter.metric=dplyr::filter(metric, ID %in% to_plot)
    
    # use for loop with rbind to filter by ID such that dataframe is ordered by sequency in which user chose IDs
    
    filter.long=dplyr::filter(long.data, policy %in% to_plot[1])
    filter.wide=dplyr::filter(wide.data, ID %in% to_plot[1])
    filter.metric=dplyr::filter(metric, ID %in% to_plot[1])
    
    
    if(!is.null(to_plot[2])){
      
      for(i in to_plot[-1]){
        
        add=dplyr::filter(long.data, policy==i)
        filter.long=rbind(filter.long, add)
        add=dplyr::filter(wide.data, ID==i)
        filter.wide=rbind(filter.wide, add)
        add=dplyr::filter(metric, ID==i)
        filter.metric=rbind(filter.metric, add)
      }
      
    }
    
    
    
    # get rank, append to data frames
    
    correction=ifelse(preferred_direction == 'min', 1, -1) # to handle metrics that should be minimized and maximizied accordingly in the ranking
    decreasing=ifelse(preferred_direction == 'min', F, T) # needed for add_lines in plot function
    
    
    if(metric_label=="order"){
      
      filter.metric$order=1:nrow(filter.metric) # add a column that simply indicates the order by which user added metrics
      # filter.metric=left_join(x=data.frame(order=as.integer(to_plot)), y=filter.metric, by="order") # reorder by the order in which user entered them
      
    }
    
    
    filter.metric$rank=rank(correction*filter.metric[[metric_label]], ties.method = 'first')
    filter.long$rank=rep(filter.metric$rank, each=length(unique(filter.long$Tier)))
    filter.wide$rank=filter.metric$rank
    
    
    ############################# plotting ################################
    text_size=labelsize
    
    n_policies=nrow(filter.wide)
    
    bar_plot=ggplot()+
      geom_bar(data=filter.long, aes(fill=Tier, x=rank, y=delta), position = 'stack', stat = 'identity', color='darkgrey')+
      geom_text(data=filter.long, aes(x=rank, y= elevation, label=v_lab),color='black', nudge_y = -4, size=text_size, check_overlap = T)+
      geom_text(data=filter.wide, aes(x=rank, y= policy_lab_y, label=policy_lab), nudge_y = 4, size=text_size, check_overlap = T)+
      # geom_text(data=filter.wide, aes(x=rank, y= policy_lab_y, label=SOM_node), nudge_y = 10, size=text_size, check_overlap = T)+ I have removed SOM node for NOW
      # scale_fill_brewer(palette='RdYlBu', direction = -1)+
      scale_fill_brewer(palette='YlOrRd', direction = 1)+
      xlab(paste(metric_label, 'rank', sep=' '))+
      ylab('pool elevation [ft msl]')+
      theme(plot.title = element_text(size=10), plot.margin=margin(t=0, r=0, b=0, l=0, unit='pt'))+
      coord_cartesian(ylim=c(895,1110), xlim=c(.5, min((n_policies+.5),20)))
    
    bar_plot=addSmallLegend(bar_plot)
    # ggtitle(paste(metric_label, 'rank for selected policies', sep=' '))+
    
    int_plot=ggplotly(p=bar_plot,height=height, tooltip = c('rank','elevation', 'volume'), dynamicTicks=T, originalData=F) # convert ggplot to interactive plotly html
    # add legend title in correct location
    int_plot=int_plot %>%   layout(legend = list(
      orientation = "v", title=list(text=" Tier "))
    )
    
    if (y_axis2==T){ # only add second y axis if y_axis2==T
      
      int_plot_2y=int_plot %>%
        add_lines(data=filter.metric, x=~sort(filter.metric$rank), y=~sort(filter.metric[[metric_label]], decreasing = decreasing), yaxis='y2',
                  inherit=FALSE, showlegend=FALSE, line=list(color='purple', width=2, dash='dash')) %>%
        layout(yaxis2 = list(overlaying = "y", side = "right",
                             tickfont = list(color = 'purple', size=10), color = 'purple',
                             title = metric_label),
               legend = list(x = 1.00, y = 1.00), xaxis=list(range=c(.5, min((n_policies+.5),20)), showticklabels=F, showline=F), yaxis=list(range=c(885,1110))
        )
      
    } else { # do not add second axis. Used on sensitivity analysis page
      
      int_plot_2y=int_plot %>%
        layout(legend = list(x = 1.00, y = 1.00), xaxis=list(range=c(.5, min((n_policies+.5),20)), title="", showticklabels=F), yaxis=list(range=c(885,1110)))
      
    }
    
    
    if (show_legend==F){
      int_plot_2y=layout(int_plot_2y, showlegend=F)
    }
    
    
    int_plot_2y$x$layout$xaxis$autorange = FALSE # need to tell plotly to NOT change the axis range to fit all data
    int_plot_2y$x$layout$yaxis$autorange = FALSE
    return(int_plot_2y)
    
    
    
  }
  
###################### custom function for ggpairs #######################
  #########################################################################
  
  # ggpairs function
  # SD_plot=function(data, mapping, dim_names=dim_names, COI=COI, box_lims=box_lims){
  #   
  #   xName=as_string(mapping$x[[2]])
  #   yName=as_string(mapping$y[[2]])
  #   
  #   ggplot(data, mapping)+
  #     geom_point(alpha=0.3, size=4.5)+
  #     annotate(geom="rect", xmin=box_lims[[xName]][1], xmax=box_lims[[xName]][2],
  #              ymin=box_lims[[yName]][1], ymax=box_lims[[yName]][2], colour="black", alpha=0 )
  #   
  #   
  # }
  
  SD_plot=function(data, mapping, dim_names=dim_names, COI=COI, box_lims=box_lims){
    
    xName=as_string(mapping$x[[2]])
    yName=as_string(mapping$y[[2]])
    
    # cat(file=stderr(), 'x,y:', c(xName, yName), "\n")
    
    xmin=box_lims[[xName]][1]
    xmax=box_lims[[xName]][2]
    
    ymin=box_lims[[yName]][1]
    ymax=box_lims[[yName]][2]
    
    my_cols=c("No"="blue", "Yes"="red")
    
    if(is.na(xmin) & is.na(ymin)){ # if no constraints, don't draw a box
      
      ggplot(data, mapping)+
        geom_point(alpha=0.3, size=3.5)+
        scale_colour_manual(values=my_cols)
      
    } else if(is.na(ymin) & !is.na(xmin)){ # constraint on x, not y
      
      ggplot(data, mapping)+
        geom_point(alpha=0.3, size=3.5)+
        geom_vline(xintercept=xmin)+
        geom_vline(xintercept=xmax)+
        scale_colour_manual(values=my_cols)
      
    } else if (!is.na(ymin) & is.na(xmin)){ # constraint on y, not x
      
      ggplot(data, mapping)+
        geom_point(alpha=0.3, size=3.5)+
        geom_hline(yintercept=ymin)+
        geom_hline(yintercept=ymax)+
        scale_colour_manual(values=my_cols)
      
    } else { # constraint on x and y, draw rectangle
      
      ggplot(data, mapping)+
        geom_point(alpha=0.3, size=3.5)+
        annotate(geom="rect", xmin=box_lims[[xName]][1], xmax=box_lims[[xName]][2],
                 ymin=box_lims[[yName]][1], ymax=box_lims[[yName]][2], colour="black", alpha=0 )+
        scale_colour_manual(values=my_cols)
      
    }
    
    
    
    
  }

  
  SD_plot2=function(data, mapping, COI=COI, box_lims=box_lims, npolicy, name_vec){
    
    
    # if (is.data.frame(box_lims)){
    #   box_lims=list(box_lims)
    # }
    
    xName=as_string(mapping$x[[2]])
    yName=as_string(mapping$y[[2]])
    
    #cat(file=stderr(), 'x,y:', c(xName, yName), "\n")
    

    if(npolicy==1){ # 1 policy
      
      box_cols=c("black", "red", "purple", "green", "brown", "orange", "yellow", "black")
      my_cols=c("No"="blue", "Yes"="red")
      my_shapes=NULL
      
    } else { # 2 policies
      
      
      box_cols=c("black", "red", "purple", "green", "brown", "orange", "yellow", "black")
      my_cols=c("First policy vulnerable"="purple", "Second policy vulnerable"="green", "Both vulnerable"= "red", "Neither vulnerable"= "blue")
      my_shapes=c("First policy vulnerable"=15, "Second policy vulnerable"=18, "Both vulnerable"= 17, "Neither vulnerable"= 20)
      
    }
    
    my_size=1
    
    # cat(file=stderr(), 'box_lims', length(box_lims), "\n")


    # draw scatter plot, no scenarios 
    p=ggplot(data, mapping)+
      geom_point(alpha=0.3, size=5)+
      scale_colour_manual(values=my_cols)+
      scale_shape_manual(values=my_shapes)+
      scale_fill_manual(values=my_cols)+
      new_scale("fill")+
      new_scale("color")
    
    for (i in 1:length(box_lims)){ # add scenarios
      xTF= xName %in% colnames(box_lims[[i]]) # true if x metric is contained in scenario
      yTF= yName %in% colnames(box_lims[[i]]) # true if y metric is contained in scenario
      
      if (xTF==F & yTF==F){ # neither x nor y is contained in this scenario
        # plot no scenarios
      } else if (xTF & yTF ==F){ # vertical line only
        
        xmin=box_lims[[i]][[xName]][1]
        xmax=box_lims[[i]][[xName]][2]
        
        p=p+
          geom_vline(xintercept=xmin, color=box_cols[i], size=my_size)+
          geom_vline(xintercept=xmax, color=box_cols[i], size=my_size)#+
          # geom_text(x=1.009*xmin, y=1.1*min(data[[yName]]),label=name_vec[i],color=box_cols[i], angle=90, inherit.aes =T, show.legend = F)
        
      } else if (xTF==F & yTF){# horizontal line only
        
        ymin=box_lims[[i]][[yName]][1]
        ymax=box_lims[[i]][[yName]][2]
        
        p=p+
          geom_hline(yintercept=ymin, color=box_cols[i], size=my_size)+
          geom_hline(yintercept=ymax, color=box_cols[i], size=my_size)#+
          # geom_text(x=1.1*min(data[[xName]]), y=.95*ymax,label=name_vec[i], color=box_cols[i], inherit.aes = T, show.legend = F)
        
        
      } else { # draw box
        
        p=p+
          annotate(geom="rect", xmin=box_lims[[i]][[xName]][1], xmax=box_lims[[i]][[xName]][2],
                   ymin=box_lims[[i]][[yName]][1], ymax=box_lims[[i]][[yName]][2], colour=box_cols[i], alpha=0, size=my_size )#+
          # geom_text(x=1.009*box_lims[[i]][[xName]][1], y=box_lims[[i]][[yName]][2], label=name_vec[i], hjust="left", vjust="top", color=box_cols[i], show.legend = F)

          # geom_rect(colour=box_cols[i],xmin=box_lims[[i]][[xName]][1], xmax=box_lims[[i]][[xName]][2],
          #          ymin=box_lims[[i]][[yName]][1], ymax=box_lims[[i]][[yName]][2], aes(fill=name_vec[i]) , alpha=0, size=my_size )+
          # scale_fill_manual(values=box_cols[i], labels=name_vec[i])+
          # new_scale("fill")+new_scale_color()

          
        
        
      }

      
    }
  
    return(p)
    
  }
  
  ############ create dummy rectangle plot to steal legend from ################
  
  scenario_legend=function(scen_names=TDcalcs$scen_names, colors=c("black", "red", "purple", "green", "brown", "orange", "yellow", "pink"),
                           size=1){
    
    xl=seq(1,10, length.out = length(scen_names))
    xh=seq(10, 10, length.out = length(scen_names)) 
    yl=seq(1,10, length.out = length(scen_names))
    yh=seq(10, 10, length.out = length(scen_names))
    
    legend.df=data.frame(xl, xh, yl, yh, Scenario=scen_names)
    
    if(length(scen_names)>=length(colors)){
      scen_names=scen_names[1:length(colors)]
    } else {
      
      start_index=length(scen_names)+1
      
      for(i in start_index:length(colors)){
        
        scen_names[i]=paste0("dummy",i)
        
        
      }
    }
  
    names(colors)=scen_names

    p=ggplot(data = legend.df, mapping=aes(fill=Scenario, color=Scenario))+
      geom_rect(xmin=xl, xmax=xh, ymin=yl, ymax=yh, alpha=0, size=size)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)
    
    p=p+theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),
                          rect = element_rect(fill = "white"), legend.title=element_text(size=16), legend.text=element_text(size=14)) # make legend background transparent
    
    p=get_legend(p) #steal legend
    p=as_ggplot(p) # convert to ggplot
    
    return(p)
    
  }
  
  
  ############ modified ggpairs diagonal plot (just changing the colors) ################
  
  modified_density = function(data, mapping, ...) {
    my_cols=c("No"="blue", "Yes"="red")
    ggally_densityDiag(data, mapping, ...) + scale_fill_manual(values=my_cols)
  }
  
  
modified_points=function(data, mapping, ...){

  my_cols=c("No"="blue", "Yes"="red")

  ggplot(data, mapping)+
    geom_point(alpha=0.3, size=3.5)+
    scale_colour_manual(values=my_cols)
  
  
}

modified_points2=function(data, mapping, ...){
  my_cols=c("First policy vulnerable"="purple", "Second policy vulnerable"="green", "Both vulnerable"= "red", "Neither vulnerable"= "blue")
  my_shapes=c("First policy vulnerable"=15, "Second policy vulnerable"=18, "Both vulnerable"= 17, "Neither vulnerable"= 20)
  
  ggplot(data, mapping)+
    geom_point(alpha=0.3, size=5)+
    scale_colour_manual(values=my_cols)+
    scale_shape_manual(values=my_shapes)+
    scale_fill_manual(values=my_cols)
  
}

modified_density2 = function(data, mapping, ...) {

  my_cols=c("First policy vulnerable"="purple", "Second policy vulnerable"="green", "Both vulnerable"= "red", "Neither vulnerable"= "blue")

  ggally_densityDiag(data, mapping, ...)+
    scale_fill_manual(values=my_cols)
}



########################### function to remove extra ggpairs labels when diagonal not wanted ##############

# taken directly from here: https://stackoverflow.com/questions/42654928/how-to-show-only-the-lower-triangle-in-ggpairs/42656454

gpairs_lower <- function(g){
  g$plots <- g$plots[-(1:g$nrow)]
  g$yAxisLabels <- g$yAxisLabels[-1]
  g$nrow <- g$nrow -1
  
  g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
  g$xAxisLabels <- g$xAxisLabels[-g$ncol]
  g$ncol <- g$ncol - 1
  
  g
}
  
  
  
  ######################### uncertainty metric range plot ########################
  
  units=read.csv("uncertainty metric units.csv", header=T)[,-1]
  
  range_plot=function(box_lims, x=x$x, units=units, annotation_size=input$RangePlotFont){
    
    scaled=matrix(nrow=nrow(box_lims), ncol=ncol(box_lims))
    for (i in 1:ncol(box_lims)){
      
      metric_name=colnames(box_lims)[i]
      
      max=max(x[[metric_name]])
      
      min=min(x[[metric_name]])
      
      scaled[,i]=(box_lims[,i]-min)/(max-min)

    }
    
    scaled[2,]=ifelse(scaled[2,]>1, 1, scaled[2,]) # if upper bound is greater than 1, set to 1
    scaled[1,]=ifelse(scaled[1,]<0, 0, scaled[1,]) # if upper bound is less than 0, set to 0
    
    scaled=data.frame(scaled)
    colnames(scaled)=colnames(box_lims)
    
    row.names(scaled)=c("lower", "upper")
    
    to_plot=data.frame(t(scaled))
    to_plot$metric=row.names(to_plot)
    to_plot$y=1:nrow(to_plot)
    
    lower_text=paste(round(unlist(box_lims[1,]),2), unlist(units[,to_plot$metric]))
    upper_text=paste(round(unlist(box_lims[2,]),2), unlist(units[,to_plot$metric]))
    
    
    p=ggplot(data=to_plot)+
      geom_linerange(mapping=aes(xmin=lower, xmax=upper, y=metric ), size=8, color="red", alpha=0.75)+
      geom_text(label=lower_text, mapping=aes(x=lower-.05, y=metric, angle=90), inherit.aes = FALSE, size=annotation_size)+
      geom_text(label=upper_text, mapping=aes(x=upper+.05, y=metric, angle=-90), inherit.aes = FALSE, size=annotation_size)+
      scale_y_discrete(name="")+
      scale_x_continuous(name="scaled metric value", limits=c(-.05,1.05))+
      ggtitle("")+
      theme(rect = element_rect(fill = "transparent")
            , panel.border = element_rect(fill="transparent"),
            plot.background = element_rect(color="transparent"), text=element_text(size=20),
            axis.text.y = element_text(angle =0, hjust = 1, debug = F, vjust=0)
            )+
      scale_y_discrete(expand = expansion(add=1)) # to make sure metric labels are not cutoff on top and bottom of plot
     
    
    return(p)
    
  }
  
  
##################### violin plot wrapper function #######################
  
  violin_plot=function(data, objective, units, source, line_y=max(data[,objective]), add_traces=F,firstpolicy
                       , policyID, objectives_all=obj_all, color_vec, addKey=FALSE, key_scale){
    
    
    if (is.null(policyID)){
      IDs=firstpolicy
    } else {
      IDs=c(firstpolicy, policyID)
    }
    
    # use for loop with rbind to filter by ID such that dataframe is ordered by sequency in which user chose IDs
    
    obj=dplyr::filter(data, ID == firstpolicy)
    
    if(!is.null(policyID)){
      
      for(i in policyID){
        
        add=dplyr::filter(obj_all, ID==i)
        
        obj=rbind(obj, add)
        
      }
      
    }
    

    # obj=dplyr::filter(obj_all, ID %in% IDs) # don't use this method because you lose the order in which user selected policies
    
    obj$ID=factor(obj$ID)
    
    
    line <- list(
      type = "line",
      line = list(color = "black", dash="dash", width=4),
      xref = "paper",
      yref = "y",
      y0=line_y,y1=line_y,
      x0=0,
      x1=1,
      layer="above"
    )
    
    
    image=list(
      source=base64enc::dataURI(file = "images/Vulnerability threshold.png"),
      layer="above",
      xref="paper",
      x=0.02,
      yref="paper",
      y=.85,
      sizex = key_scale, sizey = key_scale*max(obj[,objective]), yanchor="middle"
    )
    
    a=ggplot()+ # in geom_violin, use in_order() such that x is mapped not numerically but by order in which user selected IDs
      geom_violin(data=obj[,c("ID","TraceNumber",objective)], mapping=aes(color=ID,y=get(objective), x=fct_inorder(ID)),trim=T, fill="transparent",alpha=.4, size=.6)+
      xlab("")+
      ylab("TEST")

    ydata=seq(min(obj[,objective]),max(obj[,objective]), length.out = nrow(obj))
    
    for(i in 1:length(IDs)){
      
      if(i==1){
        line.df=data.frame(ydata, ID=IDs[i])
      } else {
        
        add.df=data.frame(ydata, ID=IDs[i])
        
        line.df=rbind(line.df, add.df)
        
      }
      
    }
    
    line.df$ID=factor(line.df$ID)
    
    a=a+
      geom_point(data=line.df, aes(x=fct_inorder(ID), y=ydata, color=ID), alpha=0)
 
    a=change_palette(a, palette = c("green", "purple", "orange", "blue", "red"))
    
    a=ggplotly(a, source = source)
    
    temp_data=dplyr::filter(obj, ID == IDs[1])
    
    a=add_trace(a, 
                type="box",
                y=temp_data[,objective],
                x=factor(1),
                width=0.2,
                line=list(color="black", width=3),
                fillcolor="transparent",
                boxmean=TRUE,
                marker=list(outliercolor="black", color="black", line=list(color="black"), opacity=0.2)
                
    )
    
    
    if(!is.null(policyID)){ # add traces
      
      for(i in 1:length(policyID)){
        
        temp_data=dplyr::filter(obj, ID == policyID[i])
        
        # add box plot
        a=add_trace(a, 
                    type="box",
                    y=temp_data[,objective],
                    x=factor(i+1),
                    width=0.2,
                    line=list(color="black", width=3),
                    fillcolor="transparent",
                    boxmean=TRUE,
                    marker=list(outliercolor="black", color="black", line=list(color="black"), opacity=0.2)
                    
        )

        
      }
      
    }
    
                
                
   
    
    if(addKey==TRUE){
      a=a %>% layout(images=image)
    }
    
    a= a %>% layout(margin=list(b=60,t=35,pad=0),
                    annotations=list(text=paste0(objective, " [",units, "]"), y=0, yref="paper",
                                     x=0, xref="paper", font=list(size=15), showarrow=F, xanchor="left"
                                     ), dynamicTicks=T
                    )
 
    a$x$layout$shapes <- line

    a=a %>% config(edits=list(shapePosition=TRUE))
    
    n=length(IDs)
    keep=(n+1):(n+n)
    remove=(1:(3*n))[-keep]

    a=style(a, hoverinfo="skip", traces = remove) # don't show hover info for violin and box plot traces

    ydata=seq(min(obj[,objective]),max(obj[,objective]), length.out = nrow(obj) )
    text_y=paste(round(ydata, 3), rep(units, nrow(obj)))

    a=style(a, text=text_y, traces = keep)
    
    return(a)
    
  }
    
  
  
  
  
  
# violin_plot=function(data, objective, units, source, line_y=max(data[,objective]), add_traces=F,firstpolicy, policyID, objectives_all=obj_all, color_vec){
#   
#   line <- list(
#     type = "line",
#     line = list(color = "black", dash="dash", width=4),
#     xref = "paper",
#     yref = "y",
#     y0=line_y,y1=line_y,
#     x0=0,
#     x1=1,
#     layer="below"
#   )
# 
#   
#   # text=list(
#   #   
#   #   text="vulnerable<br> <br>not vulnerable",
#   #   font=list(size=14),
#   #   xref = "paper",
#   #   yref = "y",
#   #   x=1,
#   #   y=line_y,
#   #   arrowside="none",
#   #   arrowidth=6,
#   #   axref="paper", ax=0,
#   #   ayref="pixel", ay=0,
#   #   xanchor="auto",
#   #   xshift=50
#   #   
#   # )
#   
#   image=list(
#     source=base64enc::dataURI(file = "images/Vulnerability threshold.png"),
#     layer="above",
#     xref="paper",
#     x=0.02,
#     yref="paper",
#     y=.85,
#     sizex = .3, sizey = 0.3*max(data[,objective]), yanchor="middle"
#   )
# 
# 
#   
#   fig<- plot_ly(
#       y = data[,objective],# x=x,
#       type = 'violin', hoverinfo="y", showlegend=F, spanmode="hard",source = source,
#       points="all", hoveron="points", pointpos=0, jitter=0, marker=list(opacity=.015, color="black"),
#       legendgroup=as.character(firstpolicy),
#       box = list(
#         visible = T
#       ),
#       meanline = list(
#         visible = T
#       ),
#       x0 = paste0(objective, " [",units, "]"),
#       color=I("green"),
#       offsetgroup="1",
#       alignmentgroup="1",
#       xaxis="x1",
#       yaxis="y1"
#     ) 
#   
#   
#   if(add_traces==T){
# 
#     for(i in 1:length(policyID)){
# 
#       if(policyID[i]=="all"){
#         IDs=1:463
#       } else {
#         IDs=policyID[i]
#       }
#       
#       alt_data=dplyr::filter(objectives_all, ID %in% IDs)
#       col=color_vec[i]
# 
# 
#       fig=fig %>% add_trace(
#         type = 'violin',
#         y = alt_data[,objective],
#         name=as.character(policyID[i]),
#         legendgroup=as.character(policyID[i]),
#         box = list(
#           visible = T
#         ),
#         meanline = list(
#           visible = T
#         ),
# 
#         color=I(col),
#         offsetgroup=as.character(i+1),
#         alignmentgroup="1",
#         xaxis="x1",
#         yaxis="y1"
#       )
# 
# 
# 
#     }
# 
#   }
# 
#   if(add_traces==T){
#     mode="group"
#   } else {
#     mode="overlay"
#   }
#   
#   
#   fig <- fig %>%
#     layout(
#       yaxis = list(
#         title = "",
#         zeroline = F
#       ), margin=list(b=50,pad=0), shapes=line,
#       violinmode = mode
#       )
#   
#   fig= fig %>% layout(
#     images=image
#   )
#   
#   fig=fig %>% config(edits=list(shapePosition=TRUE))
#   return(fig)
#   
#   
# }
  

########## violin plot subplots ############################
  
  names=colnames(obj_all)[-c(1,2)]
  my_units=c("KAF", "KAF", "%", "%", "%", "Yrs", "%", "KAF", "MAF")
  
violin_subplots=function(obj=obj$filter ,obj_names=names, unit_vec=my_units, source="VP", y_vec=thresholds$thresholds, policyID=input$policyIDcompare, color_vec,
                         add_traces=F, objectives_all=obj_all, firstpolicy, img_height=input$height, img_width=input$width, img_scale=input$scale,
                         img_format=input$format){
  
  violin_list=list()
  
  for (i in 1:length(obj_names)){
    
    
    if (i==1){ # add legend to first plot only
      
      if (length(obj_names)>1){key=0.45} else {key=0.3} # set size of key depending on how many objectives are plotted
      
      violin_list[[i]]=violin_plot(data=obj, objective = obj_names[i],units = unit_vec[i], source=source, line_y = y_vec[i], policyID=policyID, color_vec=color_vec, objectives_all=obj_all,
                                   add_traces = add_traces, firstpolicy = firstpolicy, addKey=TRUE, key_scale=key) 
    } else {
      
      violin_list[[i]]=violin_plot(data=obj, objective = obj_names[i],units = unit_vec[i], source=source, line_y = y_vec[i], policyID=policyID, color_vec=color_vec, objectives_all=obj_all,
                                   add_traces = add_traces, firstpolicy = firstpolicy, key_scale=0.3)
    }
    
    
    
    
  }
  
  n=length(obj_names)
  nrow=if(n %in% 1:2 ){1} else if (n %in% 3:4){2} else {3}

  sub=subplot(violin_list, nrows = nrow, shareX = F, titleY=F, titleX = T)
  
  # see these sources for wayt to change the download image behavior. You can change height, width, and fyle type to obtain MUCH IMPROVED quality
  # https://plotly.com/python/configuration-options/
  # https://www.rdocumentation.org/packages/plotly/versions/4.9.3/topics/config
  # https://github.com/plotly/plotly.js/blob/master/src/plot_api/plot_config.js
  
  sub=config(sub, toImageButtonOptions=list(format=img_format, height=img_height, width=img_width, scale=img_scale))
  
  
  return(hide_legend(sub))
  
}
