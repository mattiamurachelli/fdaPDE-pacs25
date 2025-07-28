# Switzerland Rainfall - Bernardi (2018)

library(fdaPDE)
rm(list=ls())
#### Load the data and build the objects to perform the smoothing ############

scale_domain <- 20000
scale_data <- 200

# load the observations 
# the 4th column contains the observations "z_i"
# the 2nd and 3rd columns contain the coordinates of the observations "p_i"
data_obs <- read.table('dati_SIC97_full.txt', header=FALSE, sep=',')
# load the coordinates of the boundary of the domain
boundary <- read.table('extended_boundary.txt', header=TRUE)

data <- as.numeric(data_obs[,4])/scale_data
range(data)

data_locations <- matrix(NA,nrow=dim(data_obs)[1],ncol=2)
data_locations[,1] <- as.numeric(data_obs[,2])/scale_domain
data_locations[,2] <- as.numeric(data_obs[,3])/scale_domain

rm(data_obs)

# create the matrix p with all the data locations and boundary locations
p <- matrix(data=NA,nrow=length(data)+dim(boundary)[1],ncol=2)
p[,1] <- c(data_locations[,1],boundary[,1])
p[,2] <- c(data_locations[,2],boundary[,2])
plot(p,pch=16,cex=0.7,col=c(rep('red',dim(data_locations)[1]),rep('black',dim(boundary)[1])))

# create the boundary segments (each row corresponds to an edge, that goes from
# the vertex in the first column to the one in the second column)
isboundary <- matrix(data=NA,nrow=dim(boundary)[1],ncol=2)
isboundary[,1] <- (length(data)+1):(length(data)+dim(boundary)[1])
isboundary[,2] <- c((length(data)+2):(length(data)+dim(boundary)[1]),(length(data)+1))

# create and plot the mesh
mesh_1 <- create.mesh.2D(p, order = 1, segments = isboundary)
mesh <- refine.mesh.2D(mesh_1, maximum_area=0.35, delaunay=FALSE)

plot(mesh,asp=1, pch=".")
box()
points(mesh$nodes[which(mesh$nodesmarkers==0),], pch=16,cex=1.2)
points(mesh$nodes[which(mesh$nodesmarkers==1),], pch=16, col='red',cex=1.2)

# write.csv(format(mesh$nodes, digits=16), "points.csv")
# write.csv(format(mesh$triangles, digits=16), "cells.csv")
# write.csv(format(mesh$nodesmarkers, digits=16), "boundary.csv")

n = length(data)

plot(mesh,asp=1, pch=".")
points(data_locations, pch=16, col="red")

# write.csv(format(data_locations, digits=16), "locs.csv")
# write.csv(format(as.matrix(data), digits=16), "response.csv")

Tri <- mesh$triangles

basisobj <- create.FEM.basis(mesh)

# ------------------------------------------------------------------------------

# Compute the area of the domain (needed for the reparametrization from lambda to rho)
omega <- 0
for(i in 1:dim(Tri)[1]){
  xA <- mesh$nodes[Tri[i,1],1]
  yA <- mesh$nodes[Tri[i,1],2]
  xB <- mesh$nodes[Tri[i,2],1]
  yB <- mesh$nodes[Tri[i,2],2]
  xC <- mesh$nodes[Tri[i,3],1]
  yC <- mesh$nodes[Tri[i,3],2]
  area_triangolo <- abs(xA*(yB-yC)+xB*(yC-yA)+xC*(yA-yB))/2
  omega <- omega + area_triangolo
}
rm(i)

# rho sequence for the GCV
rhoseq <- 10^seq(-7,-0.01,by=0.1) # 10^seq(-5,-2,by=0.05) # 10^seq(-11,-1,by=0.5)

# rho sequence for iterative optimization
rhovec <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
rhovec_l <- length(rhovec)

# graphical parameters for the plots
zlim <- range(data)
levels <- seq(from=min(data),to=max(data),length=10)

# plot data
nodes <- data_locations
z <- data

colore <- heat.colors(130)[2:101]
spessore <- seq(from=0.5,to=2.5,length=100)

Min <- min(z)
Max <- max(z)

# pdf(file = "data.pdf")
# plot(boundary,col='white',asp=1,xlab='',ylab='', xaxt="n", yaxt="n",
#      frame.plot=FALSE)
# points(nodes[which(round((z-Min)/(Max-Min)*100)==0),1],nodes[which(round((z-Min)/(Max-Min)*100)==0),2],col=heat.colors(130)[1],
#        cex=0.5-diff(seq(from=0.5,to=2,length=100))[1],pch=16)
# for(pippo in 1:100)
#   points(nodes[which(round((z-Min)/(Max-Min)*100)==pippo),1],nodes[which(round((z-Min)/(Max-Min)*100)==pippo),2],col=colore[pippo],cex=spessore[pippo],pch=16)
# 
# points(boundary,type = 'l',lwd=2)
# points(boundary[c(dim(boundary)[1],1),],type = 'l',lwd=2)
# dev.off()


# array to save the results
result_optim_LBFGSB <- array(NA,dim=c(rhovec_l,2)) # angolo e intensit?

GCV_optim <- rep(NA,rhovec_l)
GCV_smooth <- rep(NA,rhovec_l)

rho_smooth <- rep(NA,rhovec_l)

flag_rho_smooth_extreme_grid <- rep(0,rhovec_l)

b = array(0, c(2, nrow(points)))
c = 0
# K has to be estimated therefore is defined later

function_to_optimize = function(x){
  
  # x is the vector to be estimated, contains the direction of anisotropy (angle)
  # and the ratio between the two eigenvalues (intensity)
  angle <- x[1]
  intensity <- x[2]
  
  Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
  sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
  kappa <- Q%*%sigma%*%solve(Q)
  
  PDE_parameters <- list(K = kappa, b = b, c = c)
  
  lambdapar =  rho/(1-rho)*n/omega
  #print(paste0("lambda ", lambdapar))
  invisible(capture.output(suppressWarnings(suppressMessages(smoothing <- smooth.FEM(locations=data_locations, observations=data, 
                          FEMbasis=basisobj, lambda= rho/(1-rho)*n/omega, 
                          PDE_parameters = PDE_parameters,
                          DOF.evaluation=NULL)))))
  
  
  H <- sum(( eval.FEM(smoothing$fit.FEM, data_locations) - data )^2)
  #print(paste0("f(x) = ", H))
  return(H)
}  

xx <- c(pi/2,5.)
angle <- xx[1]
intensity <- xx[2]

Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
kappa <- Q%*%sigma%*%solve(Q)

print(Q)
print(sigma)
print(kappa)

start_time <- Sys.time()

# iterative optimization for each value of rho in rhovec (span (0,1), starting 
# from 0.01 and incrising to 0.99)
for(i in 1:rhovec_l){
  
  rho <- rhovec[i]
  
  if(i==1) start <- c(pi/2,5) # initial value of the parameter x
  if(i!=1) start <- result_optim_LBFGSB[i-1,] # for subsequent iteration start from the previous estimated value
  
  print(paste0("start ", start[1], " ", start[2]))
  # optim is the function that do the optimization in x
  # takes in input: - the starting value for the optimization (par)
  #                 - the function to optimize (fn)
  #                 - other parameters that regards the optimization procedure (see help("optim"))
  result_optim_LBFGSB[i,] <- optim(par = start, fn = function_to_optimize, method = "L-BFGS-B", 
                                   lower = c(0,1e-3), upper = c(pi,1000))$par
  
  # deal with the fact that the angle is "circular", restart the optimization if it ends at a boundary of [0,pi]
  if(result_optim_LBFGSB[i,1]==pi){
    print("sfigato1")
    result_optim_LBFGSB[i,] <- optim(par = c(0 ,result_optim_LBFGSB[i,2]), 
                                     fn = function_to_optimize, method = "L-BFGS-B", 
                                     lower = c(0,1), 
                                     upper = c(pi,1000))$par
  }
  if(result_optim_LBFGSB[i,1]==0) {
    print("sfigato2")
    result_optim_LBFGSB[i,] <- optim(par = c(pi,result_optim_LBFGSB[i,2]), 
                                     fn = function_to_optimize, 
                                     method = "L-BFGS-B", 
                                     lower = c(0,1), 
                                     upper = c(pi,1000))$par
  }
  print(paste0("optim:", result_optim_LBFGSB[i,1], " ", result_optim_LBFGSB[i,2]))
  ## GCV
  angle <- result_optim_LBFGSB[i,1]
  intensity <- result_optim_LBFGSB[i,2]
  
  Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
  sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
  kappa <- Q%*%sigma%*%solve(Q)
  
  PDE_parameters <- list(K = kappa, b = b, c = c)
  
  GCV_optim[i] <- smooth.FEM(locations=data_locations, observations=data, 
                             FEMbasis=basisobj, lambda= rho/(1-rho)*n/omega, 
                             PDE_parameters = PDE_parameters, 
                             lambda.selection.criterion='grid', 
                             lambda.selection.lossfunction='GCV',
                             DOF.evaluation='exact')$optimization$GCV
  
  # For the value of the estimated parameter, perform the smoothing 
  # estimating the best rho, and save the value of the corresponding GCV in GCV_smooth[i]
  smoothing_aniso <- smooth.FEM(locations=data_locations, observations=data, 
                                FEMbasis=basisobj, 
                                PDE_parameters = PDE_parameters, 
                                lambda.selection.criterion='newton', 
                                lambda.selection.lossfunction='GCV',
                                DOF.evaluation='exact',lambda.optimization.tolerance = 0.001)
  
  GCVseq <- smoothing_aniso$optimization$GCV_vector
  plot(GCVseq)
  points(which.min(GCVseq),min(GCVseq),col='red',pch=16)
  rho_smooth[i] <- rhoseq[which.min(GCVseq)]
  if(rho_smooth[i]==min(rhoseq)||rho_smooth[i]==max(rhoseq)) flag_rho_smooth_extreme_grid[i] <- 1
  
  # # final smoothing - not necessary
  # smoothing_aniso <- smooth.FEM.PDE.sv.basis(locations=data_locations, observations=data, FEMbasis=basisobj, lambda= rho_smooth[i]/(1-rho_smooth[i])*n/omega, PDE_parameters = PDE_parameters, GCV=TRUE)
  
  GCV_smooth[i] <- smoothing_aniso$optimization$GCV
  
}
end_time <- Sys.time()
print(paste0("Tot time: ", end_time-start_time, " s"))

rm(i)

# Choose as parameter the one that gives the minimum value of GCV_smooth
index <- which.min(GCV_smooth)

angle <- result_optim_LBFGSB[index,1]
intensity <- result_optim_LBFGSB[index,2]

Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
kappa <- Q%*%sigma%*%solve(Q)

PDE_parameters <- list(K = kappa, b = b, c = c)

# estimate the final solution
smoothing_final_aniso <- smooth.FEM(locations=data_locations, observations=data, 
                                    FEMbasis=basisobj, 
                                    lambda= rho_smooth[index]/(1-rho_smooth[index])*n/omega, 
                                    PDE_parameters = PDE_parameters)

# save.image("workspace.RData")

###############################################################################################################################################
### Plot the estimated tensor and the estimated field

# load("workspace.RData")
# library(fdaPDE)


library(car) # needed only for the visualization of the ellipse

rho_smooth[index]

angle 
intensity 

# ellipse

library(Matrix)
kappa_cpp = as.matrix(readMM("kappa.mtx"))

png(file = "ellipse_bernardi_2018.png")
plot(boundary,col='white',asp=1,xlab='',ylab='',xaxt="n",yaxt="n",
     frame.plot=F)
points(nodes[which(round((z-Min)/(Max-Min)*100)==0),1],nodes[which(round((z-Min)/(Max-Min)*100)==0),2],col=heat.colors(130)[1],cex=0.5-diff(seq(from=0.5,to=2,length=100))[1],pch=16)
for(pippo in 1:100)
  points(nodes[which(round((z-Min)/(Max-Min)*100)==pippo),1],nodes[which(round((z-Min)/(Max-Min)*100)==pippo),2],col=colore[pippo],cex=spessore[pippo],pch=16)

points(boundary,type = 'l',lwd=2)
points(boundary[c(dim(boundary)[1],1),],type = 'l',lwd=2)
legend(x = "topleft", bty = "n", cex = 2,
       fill=c("blue", "green3"), legend = c("Bernardi", "fdaPDE 2.0"))
ellipse(c(0,0), kappa, 3, center.cex=0, add=TRUE, col='blue', lwd=5, grid=FALSE, asp=1) #, log="", center.pch=19, segments=51, draw=TRUE, xlab="", ylab="", fill=FALSE, fill.alpha=0.3)
ellipse(c(0,0), kappa_cpp, 3, center.cex=0, add=TRUE, col='green3', lwd=5, 
        grid=FALSE, asp=1, lty=2) #, log="", center.pch=19, segments=51, draw=TRUE, xlab="", ylab="", fill=FALSE, fill.alpha=0.3)

dev.off()

