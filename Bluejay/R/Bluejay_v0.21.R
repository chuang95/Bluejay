#' Bluejay
#'
#' R package use for classify ovarian cancer patients drug response by CA125 profile
#'
#' Required packages: ggplot2, grid
#'
#' @docType package
#' @name DrugResponse
NULL
 
#########################################
##           Drug Response             ##
#########################################
#' Classify OV patients drug response by CA125 profile
#'
#' Plot CA125 profile and drug response label
#'
#' This function classifies drug response by using CA125 history profile
#'
#' @param inputfile     	CA125 pd file
#'
#' @param patient_name     	name of the patient
#'  
#' @param Days_plot			days after the first CA125 test
#'
#' @param Day_limit			time period, this function will classify patient base on the CA125 history of these days + 21 days (3 weeks) after surgery		 
#'
#' @param CA125_limit		CA125 upper limit when plotting
#'
#' @param CA125_bound		CA125 value lower bound, we will consider this patient temporary cured if her CA125 value lower than this number
#'
#' @param dfs     			default 21 days (3 weeks) after surgery
#' 
#' @return a drug response label: non-determined, sensitive, resistant
#'
#' @export
 
Bluejay.classify <- function(inputfile,patient_name,Days_plot=1000,Day_limit=200,CA125_limit=500,CA125_bound=35,dfs=21)
{
    #read input file
    #fname=strsplit(inputfile, split='.', fixed=TRUE)
    #pname=strsplit(fname[[1]][1], split=' ', fixed=TRUE)
    Pinf =read.csv(inputfile,head=TRUE)
 
    #surgery date
    s_d = as.numeric(Pinf$s_d[!is.na(Pinf$s_d)])
    s_l = length(s_d)
    Pinf=Pinf[which(Pinf[,1]<(s_d[s_l]+ Day_limit)),]
    n_nullind=which(Pinf$drug!= "NULL")
    nnl=length(n_nullind)
    last_day=Pinf$days[n_nullind[nnl]]
    last_d=Pinf$drug[n_nullind[nnl]] 
    pd = dim(Pinf)
    
    #get therapies name and working time range
    Therapy_ds = as.character(Pinf$drug[!is.na(Pinf$drug)])
    Therapy =Therapy_ds[which(Therapy_ds!="NULL")]
    Therapy =Therapy[which(Therapy!="continue")]
 
    l=length(Therapy)
    din=list()
    for (i in 1:l)
    {
        temp=which(Pinf$drug== Therapy[i])
        din=cbind(din,t(Pinf$days[temp]))
    }
 
    din=as.numeric(din)
 
    din=sort(unique(din))
 
    if(din[l]!=last_day)
    {
        din=c(din,last_day)
        l=l+1
    }
 
    if(din[l]==last_day && last_d!="continue")
    {
        din=c(din,Pinf$days[(n_nullind[nnl]+1)])
        l=l+1
    }
 
    nd=Pinf$days[which(Pinf$drug== "NULL")]
 
    ndl=length(nd)
 
    k=1
    xstart=list()
    xend=list()
    kflag=FALSE
    for (j in 1:(l-1))
    {
 
        while(din[j]>nd[k] && !kflag)
        {
            k=k+1
            if(k>ndl)
            {
                k=ndl
                kflag=TRUE
                break
            }
        }
        if(din[j]<nd[k] && !kflag && j<(l-1))
        {
            xstart=c(xstart,din[j])
            if(din[j+1]<nd[k])
            {
                xend=c(xend,din[j+1]-1)
            }
            else
            {
                xend=c(xend,nd[k])
                k=k+1
            }
        }
        else if(kflag && din[j]>nd[k] && j<(l-1))
        {
            xstart=c(xstart,din[j])
            xend=c(xend,din[j+1]-1)
        }
        else if(j==(l-1))
        {
            xstart=c(xstart,din[j])
        }
        else
        {
            xstart=c(xstart,din[j])
            k=k+1
        }
 
        if(k>ndl)
        {
            k=ndl
            kflag=TRUE
        }
    }
    xend=c(xend,din[j+1])
    xstart=as.numeric(xstart)-(s_d[1]+dfs)
    xend=as.numeric(xend)-(s_d[1]+dfs)
    xstart[which(xstart<0)]=0
    xend[which(xend<0)]=0
    rects <- data.frame(xstart, xend, Therapy)
 
    #remove NA and plot CA125 profile
    Pinf=Pinf[!is.na(Pinf$ca125),]
    #y_l=max(Pinf[,2]+300)
    l_t=length(Therapy)
    Therapy_s= Therapy
    Therapy_carb_taxol=NA
    tct=1
    for(t in 1:l_t)
    {
        Therapy_s[t]=substr(Therapy[t],1,2)
        if(grepl("Carboplatin", Therapy[t]) || grepl("Taxol", Therapy[t]))
        {
            Therapy_carb_taxol[tct]=t
            tct=tct+1
        }
    }
 
    Pinf[,1]=Pinf[,1]-(s_d[1]+dfs)
    Pinf=Pinf[which(Pinf[,1]>0),]
 
    solid_p_x=Pinf[which(Pinf[,1]<Day_limit),]
    solid_p=Pinf[which(Pinf[,2]<= CA125_limit),]
    #Pinf[which(Pinf[,2]>CA125_limit),2]= CA125_limit
 
    s_n=length(solid_p_x[,2])
 	
 	#classify
 	label=NULL
	if(length(s_d)>1 & (s_d[2]-s_d[1])<Day_limit)
	{
		label="non-determined"
	}else
	{
		tl=length(Pinf[,1])
		if(tl<2)
		{
			label="non-determined"
		}else
		{
			average_drop_rate=Bluejay.mav(Pinf[,2],2)#movingAverage(Pinf[,2], n=s_n, centered=FALSE)
			adr=mean(diff(Bluejay.mav(solid_p_x[,2],2))/Bluejay.mav(solid_p_x[,2],2)[c(1:(s_n-1))],na.rm=TRUE)
			if(Pinf[1,2]<= CA125_bound & length(which(solid_p_x[,2]> (CA125_bound+5)))<1)
			{
				label="non-determined"
			}else if(Pinf[1,2]>= CA125_bound & (adr < (-0.3) || length(which(solid_p_x[,2] < CA125_bound))>0) )
			{
				label="sensitive"
			}else if(Pinf[1,2]>= CA125_bound & adr > (-0.3) )
			{
				label="resistant"
			}else
			{
				label="non-determined"
			}
		}
	}
 
    p=Bluejay.plot(inputfile, patient_name,label,Days_plot,CA125_limit)
    
    pdf(paste("patient_", patient_name, "_ca125_plot.pdf", sep = ""))
    print(p)
    dev.off()
    
 	return(label)
}
 
Bluejay.mav <- function(x,n=5)
{
	filter(x,rep(1/n,n), sides=2)
}

Bluejay.plot <- function(inputfile,patient_name,label,Days_plot,CA125_limit)
{
    #fname=strsplit(inputfile, split='.', fixed=TRUE)
    #pname=strsplit(fname[[1]][1], split=' ', fixed=TRUE)
    Pinf =read.csv(inputfile,head=TRUE)
	Pinf=Pinf[which(Pinf[,1]<Days_plot),]
	
	pd = dim(Pinf)

	n_nullind=which(Pinf$drug!= "NULL")
	nnl=length(n_nullind)
	last_day=Pinf$days[n_nullind[nnl]]
	last_d=Pinf$drug[n_nullind[nnl]]
	s_d = as.numeric(Pinf$s_d[!is.na(Pinf$s_d)])

	#rects <- data.frame(xstart = c(1,10,1324), xend = c(9,1323,1399), col = as.character(unique(Pinf$drug)))
	Therapy_ds = as.character(Pinf$drug[!is.na(Pinf$drug)])
	Therapy =Therapy_ds[which(Therapy_ds!="NULL")]
	Therapy =Therapy[which(Therapy!="continue")]

	l=length(Therapy)
	din=list()
	for (i in 1:l)
	{
	    temp=which(Pinf$drug== Therapy[i])
	    din=cbind(din,t(Pinf$days[temp]))
	}

	din=as.numeric(din)

	din=sort(unique(din))

	if(din[l]!=last_day)
	{
 	   din=c(din,last_day)
 	   l=l+1
	}

	if(din[l]==last_day && last_d!="continue")
	{
 	   din=c(din,Pinf$days[(n_nullind[nnl]+1)])
	    l=l+1
	}

	nd=Pinf$days[which(Pinf$drug== "NULL")]

	ndl=length(nd)

	k=1
	xstart=list()
	xend=list()
	kflag=FALSE
	for (j in 1:(l-1))
	{
    
 	   while(din[j]>nd[k] && !kflag)
 	   {
	        k=k+1
	        if(k>ndl)
	        {
	            k=ndl
	            kflag=TRUE
	            break
	        }
	    }
 	   if(din[j]<nd[k] && !kflag && j<(l-1))
 	   {
			xstart=c(xstart,din[j])
			if(din[j+1]<nd[k])
			{
				xend=c(xend,din[j+1]-1)
			}
			else
			{
	            xend=c(xend,nd[k])
	            k=k+1
	        }
	    }
	    else if(kflag && din[j]>nd[k] && j<(l-1))
	    {
	        xstart=c(xstart,din[j])
	        xend=c(xend,din[j+1]-1)
	   }
	   else if(j==(l-1))
	    {
	        xstart=c(xstart,din[j])
	    }
	    else
	    {
	        xstart=c(xstart,din[j])
	        k=k+1
	    }
    
	    if(k>ndl)
	    {
	        k=ndl
	        kflag=TRUE
	    }
	}
	xend=c(xend,din[j+1])

	xstart=as.numeric(xstart)
	xend=as.numeric(xend)

	rects <- data.frame(xstart, xend, Therapy)

	if(length(Pinf[!is.na(Pinf$ca125),2])==0)
	{
		print("No CA125 value with in Days_plot, please increase number of Days_plot")
		return(0)
	}
	
	Pinf=Pinf[!is.na(Pinf$ca125),]
	y_l=max(Pinf[,2]+300)
	l_t=length(Therapy)
	Therapy_s= Therapy

	for(t in 1:l_t)
	{
		Therapy_s[t]=substr(Therapy[t],1,2)
	}


	#solid_in=which(Pinf[,2]<=CA125_limit)
	Pinf[which(Pinf[,2]>CA125_limit),2]=CA125_limit
	solid_p=Pinf[which(Pinf[,2]<CA125_limit),]

	#pdf(paste("patient_", patient_name, "_ca125_plot.pdf", sep = ""))

	if(dim(solid_p)[1]>0){
		p <- ggplot() +geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = Therapy), alpha = 0.4)+
	geom_line(data= Pinf, aes(x=days, y=ca125, group=1)) + geom_point(data= Pinf, aes(x=days, y=ca125, group=1,shape = factor(1)),na.rm = TRUE,colour="grey50", size = 3)+guides(shape=FALSE)+ scale_shape(solid = FALSE)+geom_point(data= solid_p, 	aes(x=days, y=ca125, group=1),na.rm = TRUE,size = 3)+theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = 		element_blank())
	}else{
		p <- ggplot() +geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = Therapy), alpha = 0.4)+
	geom_line(data= Pinf, aes(x=days, y=ca125, group=1)) + geom_point(data= Pinf, aes(x=days, y=ca125, group=1,shape = factor(1)),na.rm = TRUE,colour="grey50", size = 3)+guides(shape=FALSE)+ scale_shape(solid = FALSE)+theme_bw() + 	theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
	}

	p <-p + geom_vline(xintercept = s_d,colour="red", linetype = "dotted",size = 1)+ ggtitle(paste("Patient ", patient_name, " CA125", "\n", label, sep = "")) + theme(plot.title = element_text(size=22, face="bold.italic")) +theme(legend.position = 	"top",legend.direction="vertical")+ xlim(0, Days_plot)+ylim(0,CA125_limit)

	#dev.off()
	return(p)

}
