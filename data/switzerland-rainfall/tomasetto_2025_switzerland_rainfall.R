# Switzerland Rainfall - Tomasetto (2024)
#detach("package:fdaPDE", unload = TRUE)
library(fdaPDEtomasetto) #
library(viridis)
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

# boundary anticlockwise
# plot(boundary)
# points(boundary[1,1],boundary[1,2], col="red", pch=16)
# points(boundary[2,1],boundary[2,2], col="red", pch=16)

mesh_tmp <- create.mesh.2D(boundary, segments = cbind(1:nrow(boundary), 
                                                   c(2:nrow(boundary),1)))
plot(mesh_tmp)

mesh <- refine.mesh.2D(mesh_tmp, minimum_angle = 20, maximum_area = 0.1)
dim(mesh$nodes)
plot(mesh)

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
#mesh_1 <- create.mesh.2D(p, order = 1, segments = isboundary)
#mesh <- refine.mesh.2D(mesh_1, maximum_area=0.25, minimum_angle = 30)

par(mar=c(0,0,0,0))
plot(mesh,asp=1, pch=".")
box()
points(mesh$nodes[which(mesh$nodesmarkers==0),], pch=16,cex=1.2)
points(mesh$nodes[which(mesh$nodesmarkers==1),], pch=16, col='red',cex=1.2)

Tri <- mesh$triangles

basisobj <- create.FEM.basis(mesh)

n <- length(data)

# graphical parameters for the plots
zlim <- range(data)
levels <- seq(from=min(data),to=max(data),length=10)

# plot data
nodes <- data_locations
z <- data

width <- seq(from=0.5,to=2.5,length=100)

Min <- min(z)
Max <- max(z)

color=inferno(100)
plot(boundary,col='white',asp=1,xlab='',ylab='',axes=FALSE)
points(nodes[which(round((z-Min)/(Max-Min)*100)==0),1],nodes[which(round((z-Min)/(Max-Min)*100)==0),2],col=heat.colors(130)[1],cex=0.5-diff(seq(from=0.5,to=2,length=100))[1],pch=16)
for(i in 1:100)
  points(nodes[which(round((z-Min)/(Max-Min)*100)==i),1],nodes[which(round((z-Min)/(Max-Min)*100)==i),2],col=color[i],cex=width[i],pch=16)

points(boundary,type = 'l',lwd=2)
points(boundary[c(dim(boundary)[1],1),],type = 'l',lwd=2)

###############################################################################
### anisotropic smoothing
###############################################################################
b = array(0, c(2, nrow(points)))
c = 0
angle <- pi/2
intensity <- 5
Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
kappa <- Q%*%sigma%*%solve(Q)
parameter_cascading = list(diffusion = c('K','L-BFGS-B'), advection = c('b','L-BFGS-B'))

PDE_parameters <- list(K=kappa, b=b, c=c, parameter_cascading = parameter_cascading)

start_time <- Sys.time()
smoothing_aniso <- fdaPDEtomasetto::smooth.FEM(locations=data_locations, observations=data, 
                              FEMbasis=basisobj, PDE_parameters = PDE_parameters,
                              lambda.selection.criterion='newton', 
                              lambda.selection.lossfunction='GCV',
                              DOF.evaluation='exact')

end_time <- Sys.time()
final_time = end_time - start_time

smoothing_aniso$parameter_cascading$K_direction
smoothing_aniso$parameter_cascading$K_eigenval_ratio

smoothing_aniso$parameter_cascading$b_direction
smoothing_aniso$parameter_cascading$b_intensity

# plot ellipses ----------------------------------------------------------------
library(car) # needed only for the visualization of the ellipse

angle <- smoothing_aniso$parameter_cascading$K_direction
angle
intensity <- smoothing_aniso$parameter_cascading$K_eigenval_ratio
intensity
K <- smoothing_aniso$parameter_cascading$K # K <- cbind(c(1.22961, 1.001),c(1.001,1.6281))
b <- smoothing_aniso$parameter_cascading$b # b <- c(4.2195, 5.5072)

smoothing_aniso$parameter_cascading$b_direction # diciamo coerenti
smoothing_aniso$parameter_cascading$b_intensity # diciamo coerenti

# capiamo se almeno le simulazioni sono coerenti ...


# plot ellipses
png("ellipse_tomasetto_impl.png")
#par(bg = "transparent")
color=inferno(100)
plot(boundary,col='transparent',asp=1,xlab='',ylab='',xaxt="n",yaxt="n",axes=FALSE)
points(nodes[which(round((z-Min)/(Max-Min)*100)==0),1],nodes[which(round((z-Min)/(Max-Min)*100)==0),2],col=heat.colors(130)[1],cex=0.5-diff(seq(from=0.5,to=2,length=100))[1],pch=16)
for(i in 1:100)
  points(nodes[which(round((z-Min)/(Max-Min)*100)==i),1],nodes[which(round((z-Min)/(Max-Min)*100)==i),2],col=color[i],cex=width[i],pch=16)
points(boundary,type = 'l',lwd=2)
points(boundary[c(dim(boundary)[1],1),],type = 'l',lwd=2)
ellipse(c(0,0), K, 1, center.cex=0, add=TRUE, col='#21a7f5', lwd=6, grid=FALSE, asp=1)
arrows(x0=0,y0=0,x1=b[1],y1=b[2], col='#21a7f5',lwd=6, angle = 35, length = 0.45)
dev.off()