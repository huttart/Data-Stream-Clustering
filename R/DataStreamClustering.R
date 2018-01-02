denstream <- function(r = 0.1,mu = 1,init_points = 500,lambda = 0.001,beta = 0.2,offline = 2,filename,N,SD,MEAN,dimension = c(1,2),char_split = " ",return_cluster = FALSE)
{
	# par(pty="s")
	# plot(-5:5,seq(-5,5,length=11),type="n",xlab="",ylab="",main="Test draw.circle")
	plot(-2:2,seq(-2,2,length=5),type="n",xlab="",ylab="")

	# Denstream_model = Denstream$new(r,mu,init_points,lambda,beta,offline,N,SD,MEAN,dimension,char_split,return_cluster) 
	Denstream_model <- new(Denstream,r,mu,init_points,lambda,beta,offline)
	# stop()
	data = file(filename,"r")
	pb = txtProgressBar(min = 0, max = N, initial = 0,char = "#",style = 3) 
	n = 0
	while(1){
		n = n+1
		record_n = readLines(data, n = 1)
		setTxtProgressBar(pb,n)
		if ( length(record_n) == 0 || n > N) {
			n = n-1
		    break
		}
		tmp 	= 	strsplit(record_n,char_split)[[1]]
		tmp 	= 	tmp[ tmp != ""]
		input 	= 	as.numeric(tmp[dimension])
		input 	= 	((input - MEAN)/SD)	# Normalization
		# print(input)
		# points(input[1],input[2],col = "gray",cex = 0.48)
		Denstream_model$trainOnInstance(input)

	}
	close(data)
	micro = Denstream_model$end()
	# return(micro)

	data = file(filename,"r")
	pb = txtProgressBar(min = 0, max = N, initial = 0,char = "#",style = 3) 
	n = 0
	plot(-2:2,seq(-2,2,length=5),type="n",xlab="",ylab="")

	result = as.data.frame(matrix(0, ncol = 2, nrow = N))
	while(1){
		n = n+1
		record_n = readLines(data, n = 1)
		setTxtProgressBar(pb,n)
		if ( length(record_n) == 0 || n > N) {
			n = n-1
		    break
		}
		tmp 	= 	strsplit(record_n,char_split)[[1]]
		tmp 	= 	tmp[ tmp != ""]
		input 	= 	as.numeric(tmp[dimension])
		input 	= 	((input - MEAN)/SD)	# Normalization
		class 	=	tmp[length(tmp)]
		# print(input)
		col = reassign(input,micro)
		# print(col)
		# points(input[1],input[2],col = col)

		result[n,] =  c(class,col)

	}
	close(data)

	measurement(result)




	# return(result)
}

# MUDI <- function(lambda = 0.001,alpha = 0.2,filename,N,SD,MEAN,maxmax,dimension,char_split = " ",return_cluster = FALSE) {
MUDI <- function(lambda = 0.001,alpha = 0.6,filename,N,min,max,dimension = c(1,2),char_split = " ",return_cluster = FALSE,minpts,gridsize,numofgrid) {
	plot(0:100,seq(0,100,length=101),type="n",xlab="",ylab="",main="Test draw.circle")

	Mudi_model <- new(Mudi,numofgrid,lambda,alpha,gridsize,minpts)
	data = file(filename,"r")
	pb = txtProgressBar(min = 0, max = N, initial = 0,char = "#",style = 3) 
	n = 0
	while(1){
		n = n+1
		record_n = readLines(data, n = 1)
		setTxtProgressBar(pb,n)
		if ( length(record_n) == 0 || n > N) {
			n = n-1
		    break
		}
		tmp 	= 	strsplit(record_n,char_split)[[1]]
		tmp 	= 	tmp[ tmp != ""]
		input 	= 	as.numeric(tmp[dimension])
		# input 	= 	((input - MEAN)/SD)	# Normalization
		input 	=	((input - min)/(max - min)*(100 - 0)) + 0
		# print(input)
		Mudi_model$trainOnInstance(input)
		# stop()

	}
	close(data)
	center = Mudi_model$result()
	for(i in 1:length(center)){points(center[[i]][1],center[[i]][2],col = center[[i]][3]+1)}
	plot(0:100,seq(0,100,length=101),type="n",xlab="",ylab="")
	cat("Reassinging !!!\n")
	data = file(filename,"r")
	pb = txtProgressBar(min = 0, max = N, initial = 0,char = "#",style = 3) 
	n = 0
	result = as.data.frame(matrix(0, ncol = 2, nrow = N))
	while(1){
		n = n+1
		record_n = readLines(data, n = 1)
		setTxtProgressBar(pb,n)
		if ( length(record_n) == 0 || n > N) {
			n = n-1
		    break
		}
		tmp 	= 	strsplit(record_n,char_split)[[1]]
		tmp 	= 	tmp[ tmp != ""]
		input 	= 	as.numeric(tmp[dimension])
		class 	=	tmp[length(tmp)]
		input 	=	((input - min)/(max - min)*(100 - 0)) + 0
		# print(input)
		col = reassign(input,center)
		result[n,] =  c(class,col)
		# print(reassign(input,tmp))
		# points(input[1],input[2],col = col+1)
		# stop()

	}
	close(data)

	measurement(result)
	return(result)

}

HCM <- function(filename,N,min,max,dimension = c(1,2),char_split = " ",minpts,r,length = 3) {
	# plot(0:100,seq(0,100,length=101),type="n",xlab="",ylab="",main="Test draw.circle")
	par(pty="s")

	plot(-10:20,seq(-10,20,length=31),type="n",xlab="",ylab="",main="Test draw.circle")
	HCM_model <- new(HCMStream,r,minpts,length)
	data = file(filename,"r")
	pb = txtProgressBar(min = 0, max = N, initial = 0,char = "#",style = 3) 
	n = 0
	while(1){
		n = n+1
		record_n = readLines(data, n = 1)
		setTxtProgressBar(pb,n)
		if ( length(record_n) == 0 || n > N) {
			n = n-1
		    break
		}
		tmp 	= 	strsplit(record_n,char_split)[[1]]
		tmp 	= 	tmp[ tmp != ""]
		input 	= 	as.numeric(tmp[dimension])
		# input 	= 	((input - MEAN)/SD)	# Normalization
		# input 	=	((input - min)/(max - min)*(100 - 0)) + 0
		# print(input[0])
		# points(input[1],input[2],col = 'gray',cex = 0.5)

		HCM_model$trainOnInstance(input)
		# stop()

	}
	close(data)
	center = HCM_model$result()
	# print(center)

	x = c()
	y = c()
	g = c()
	col = c()
	for(i in 1:length(center)){
		if (length(center[[i]]) <= 3) {
			# points(center[[i]][1],center[[i]][2],col = center[[i]][3]+1)
			if(center[[i]][3] > minpts){
				draw.circle(center[[i]][1],center[[i]][2],r = r,border = 'red')
			} else {
				# draw.circle(center[[i]][1],center[[i]][2],r = r)
			}
		} else if (length(center[[i]]) > 3) {
			C = c(center[[i]][1],center[[i]][2])
			I = c(center[[i]][3],center[[i]][4])
			L = center[[i]][5]
			# print(C)

			vCL = ((I*L)+C) - C
			vR = c(-vCL[2],vCL[1])
			vR = vR/(sqrt(vR[1]^2 + vR[2]^2))
			vR = (C+(vR*r)) - C
			# print(vR)
			# print(vCL)
			point1 = C + (vCL+vR)
			vR = vR*(-1)
			point2 = C + (vCL+vR)
			vCL = vCL*(-1)
			point3 = C + (vCL+vR)
			vR = vR*(-1)
			point4 = C + (vCL+vR)
			# print(point2)
			polygon(c(point1[1],point2[1],point3[1],point4[1]),c(point1[2],point2[2],point3[2],point4[2]),border = center[[i]][7])
			# draw.circle
			x = c(x,point1[1],point2[1],point3[1],point4[1])
			y = c(y,point1[2],point2[2],point3[2],point4[2])
			g = c(g,i,i,i,i)
			col = c(col,center[[i]][7],center[[i]][7],center[[i]][7],center[[i]][7])


		}
	}

	data <- data.frame(y=y, x=x,g=g)
	# ggplot(data) + geom_polygon(aes(x=x, y=y,group=g), colour="black",fill = NA) + geom_point(aes(x=x,y=y))
	# return(data)
	ggplot(data) + geom_polygon(aes(x=x, y=y,group=g), colour=col,fill = NA,lwd = 1.1) + geom_point(aes(x=x,y=y),lwd = 2,colour=col) +  theme_bw() +
    theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank()) 
	# return(data)

	# qplot(long, lat, data = m, group = group, fill = breaks, geom = "polygon")

}

loadModule("yada",TRUE)
loadModule("MOD_Mudi",TRUE)
loadModule("MOD_HCM",TRUE)




denstream2 <- function () {
	i = 40
	while(i <= 100) {
		print(i)
		denstream(filename = "../Coursera-R programming/kddcup_data_10_percent_corrected_classNum.txt",SD = sd2,MEAN = mean2,r = 0.35,N = 500000,mu = i,char_split = "\t",dimension = 1:34)
		i = i + 10
	}
}
