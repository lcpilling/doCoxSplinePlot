
#####################################################################################################################
#####################################################################################################################
####
####  Luke C. Pilling | L.Pilling@exeter.ac.uk
####      Epidemiology and Public Health Group, University of Exeter Medical School, U.K.
####
####  Function adapted from analyses by Ambarish Dutta whilst based at the UofE Medical School, Epidemiology group
####    In particular: Dutta et al, 2013. Uric acid measurement improves prediction of cardiovascular mortality in later life. 
####                       J Am Geriatr Soc. 2013 Mar;61(3):319-26.
####                       http://www.ncbi.nlm.nih.gov/pubmed/23496291
####
####  Plots a smoothed spline curve for the hazard ratio (y) for each value of the primary exposure (x) from a 
####      Cox proportional hazard model. 
####    
####  Version 0.150324
####
#####################################################################################################################
#####################################################################################################################
####    
####  For example:
#### 
####      # create the 'survival object'
####      surv.death <- na.omit( Surv(data$age_death , data$dead) )
####            
####      # do Cox model -- can include cofactors if desired -- primary independent variable (exposure of interest) must be pspline()
####      pham.fit <- coxph( surv.death ~ pspline(data$albumin, df=4) + data$sex )
####            
####      # use doCoxSplinePlot() to plot the smoothed curve (including CI's) for the Cox model
####      doCoxSplinePlot(x=data$albumin, 
####                      fit=pham.fit, 
####                      x.lab="Albumin (g/dL)", 
####                      title="Circulating albumin and risk of mortality")
####      
#####################################################################################################################
#####################################################################################################################
####
####  Essential inputs;
####	x		==  vector of x-values (the primary exposure in the Cox model -- used with pspline() function)
####	fit		==  CoxPH object -- fitted model
####	x.lab		==  label for x-axis, e.g. "Serum albumin (g/dL)" - whatever your exposure is
####	title		==  title for the plot
####
####  Optional inputs;
####	subtitle	== additional subtitle
####                               
#####################################################################################################################
#####################################################################################################################

doCoxSplinePlot <- function(x		=stop("Need to provide vector of x data: x=..."), 
                            fit		=stop("Need to provide CoxPH model with pspline term: fit=..."),
                            x.lab	=stop("Provide an x-axis label: x.lab=..."), 
                            title	=stop("Provide a title: title=..."), 
                            subtitle	=""
                           )  
{
    
    ## check packages are loaded
    require(pspline)
    require(survival)
    
    ## get predicted values for fitted spline
    predicted <- predict(fit , type = "terms" , se.fit = TRUE , terms = 1)
    
    ## determine axis limits
    smsp <- sm.spline(x , exp(predicted$fit))
    y.lim <- c( 0, ceiling(max(smsp$ysmth)/0.5)*0.5 )
    x.lim <- c(min(na.omit(x)), max(na.omit(x)))
    
        ## ticks require rounding to nearest sensible number (that's what the plot function does)
        x.limDist  <- (max(na.omit(x)) - min(na.omit(x)))			## 'distance' or span of x variable
        x.tickDist <- 0.005							## default value to round x-limits up/down to
        if (x.limDist >= 0.01 & x.limDist < 0.1)	x.tickDist <- 0.01	## determine whether it would be better to add ticks less frequently
        if (x.limDist >= 0.1  & x.limDist < 1)		x.tickDist <- 0.05
        if (x.limDist >= 1    & x.limDist < 10)		x.tickDist <- 0.1
        if (x.limDist >= 10   & x.limDist < 100)	x.tickDist <- 1
        if (x.limDist >= 100)				x.tickDist <- 10
    
        x.tickDist.5 <- x.tickDist * 5					## for rounding use 5x the distance between the smaller ticks?
        x.tickLim <- c( floor(  min(x)/x.tickDist.5)*x.tickDist.5, 	## round x-axis limits using the distance between ticks
                        ceiling(max(x)/x.tickDist.5)*x.tickDist.5 )
    
    ## plot -- setup plot with title and labels
    par(fig=c(0, 1, 0, 1))
    plot(0 , xlab=x.lab , ylab = "Hazard Ratio" , main = title , type = "n" , xlim=x.lim , ylim=c(-0.2,y.lim[2]), cex=1.2)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="grey98")
    mtext( subtitle, cex=1.1 )     
    
    # axis additional ticks & reference line
    axis(1, at=seq(x.tickLim[1], x.tickLim[2], x.tickDist), labels=F, tck=-0.009)
    axis(2, at=seq(y.lim[1], y.lim[2], 0.25), labels=F, tck=-0.009)
    axis(2, at=seq(y.lim[1], y.lim[2], 0.5), labels=F)
    abline( h = 1 , col = "grey60" , lty = 2 , lwd = 1.5)

    # add CoxPH results/fitted lines
    lines( sm.spline(x , exp(predicted$fit)) , col = "red" , lwd = 1.5)
    lines( sm.spline(x , exp(predicted$fit + 1.96 * predicted$se)) , col = "orange" , lty = 2 , lwd = 1.5)
    lines( sm.spline(x , exp(predicted$fit - 1.96 * predicted$se)) , col = "orange" , lty = 2 , lwd = 1.5)
   
    # plot boxplot to show distribution
    par(fig=c(0,1,0,0.3), new=TRUE)
    boxplot(x, horizontal=TRUE, axes=FALSE, pch=19, cex=0.5)

}

