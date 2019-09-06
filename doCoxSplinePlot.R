
#####################################################################################################################
#####################################################################################################################
####
####  Luke C. Pilling | L.Pilling@exeter.ac.uk
####      Epidemiology and Public Health Group, University of Exeter Medical School, U.K.
####
####  Function adapted from analyses by Ambarish Dutta 
####          Dutta et al, 2013. Uric acid measurement improves prediction of cardiovascular mortality in later life. 
####          J Am Geriatr Soc. 2013 Mar;61(3):319-26.
####          http://www.ncbi.nlm.nih.gov/pubmed/23496291
####
####  Plots a smoothed spline curve for the hazard ratio (y) for each value of the primary exposure (x) from a 
####      Cox proportional hazard model. 
####    
####  Version 0.190611
####
#####################################################################################################################
#####################################################################################################################
####    
####  For example:
#### 
####      # load survival & psline packages
####      library(survival)
####      library(pspline)
#### 
####      # create the 'survival object'
####      surv.death <- na.omit( Surv(data$age_death , data$dead) )
####            
####      # do Cox model -- including covariates -- primary independent variable (exposure of interest) must be pspline()
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
####	x       ==  vector of x-values (the primary exposure in the Cox model -- used with pspline() function)
####	fit     ==  CoxPH object -- fitted model
####	x.lab   ==  label for x-axis, e.g. "Serum albumin (g/dL)" - whatever your exposure is
####	title   ==  title for the plot
####
####  Optional inputs;
####	subtitle	    == additional subtitle
####    returnValues    == return a data.frame of the x and y values (default is FALSE)
####    legendPos       == default is NO results legend. Specify "topleft" or similar to show the results on the plot
####    distBoxPlot     == show bowplot of distribution? Default is FALSE
####    yMax            == manually set the y-axis maximum? (minimum always 0)
####    xLims           == manually set the x-axis limits? (needs to be a vector of length 2 - min or max can be NULL)
####    bg_box_col      == specific background colour of plot (default = grey98)
####    null_line_type  == horizontal reference line type (default = 2)
####    null_line_col   == horizontal reference line colour (default = grey60)
####    null_line_wd    == horizontal reference line width (default = 1.5)
####    h_lines_pos     == vector of y axis positions for grid lines if desired (default = NULL)
####    h_lines_type    == type of lines for y axis (default = 2)
####    h_lines_col     == colour of lines for y axis (default = grey60)
####    h_lines_wd      == width of lines for y axis (default = 1)
####    v_lines_pos     == vector of x axis positions for grid lines if desired (default = NULL)
####    v_lines_type    == type of lines for x axis (default = 2)
####    v_lines_col     == colour of lines for x axis (default = grey60)
####    v_lines_wd      == width of lines for x axis (default = 1)
####    ...             == other plot options (e.g. cex)
####                               
#####################################################################################################################
#####################################################################################################################

require(pspline)
require(survival)

doCoxSplinePlot <- function(x		=stop("Need to provide vector of x data: x=..."), 
                            fit		=stop("Need to provide CoxPH model with pspline term: fit=..."),
                            x.lab	=stop("Provide an x-axis label: x.lab=..."), 
                            title	=stop("Provide a title: title=..."), 
                            subtitle	="",
                            returnValues=FALSE,
                            legendPos	=NULL,
                            distBoxPlot	=FALSE,
                            yMax		=NULL,
							xLims		=c(NULL,NULL),
							bg_box_col	="grey98",
							null_line_type=2,
							null_line_col ="grey60",
							null_line_wd  =1.5,
							h_lines_pos   =NULL,
							h_lines_type  =2,
							h_lines_col   ="grey60",
							h_lines_wd    =1,
							v_lines_pos   =NULL,
							v_lines_type  =2,
							v_lines_col   ="grey60",
							v_lines_wd    =1,
                            ...
                           )  
{
    
    ## check packages are loaded
    require(pspline)
    require(survival)
    
    ## get predicted values for fitted spline
    predicted <- predict(fit , type = "terms" , se.fit = TRUE , terms = 1)
    
    ## determine axis limits
    smsp     <- sm.spline(x , exp(predicted$fit))
    smsp.u95 <- sm.spline(x , exp(predicted$fit + 1.96 * predicted$se))
    smsp.l95 <- sm.spline(x , exp(predicted$fit - 1.96 * predicted$se))
    y.lim <- c( 0, ceiling(max(smsp$ysmth)/0.5)*0.5 )
    x.lim <- c(min(na.omit(x)), max(na.omit(x)))
    
	## have limits been manually specified?
	if ( !is.null(yMax) )      y.lim[2] <- yMax
	if ( !is.null(xLims[1]) )  x.lim[1] <- xLims[1]
	if ( !is.null(xLims[2]) )  x.lim[2] <- xLims[2]
	
        ## get x ticks by rounding to nearest sensible number 
        x.limDist  <- (x.lim[2] - x.lim[1])			## 'distance' or span of x axis
        x.tickDist <- 0.005							## default value to round x-limits up/down to
        if (x.limDist >= 0.01 & x.limDist < 0.1)	x.tickDist <- 0.01	## determine whether it would be better to add ticks less frequently
        if (x.limDist >= 0.1  & x.limDist < 1)		x.tickDist <- 0.05
        if (x.limDist >= 1    & x.limDist < 10)		x.tickDist <- 0.1
        if (x.limDist >= 10   & x.limDist < 100)	x.tickDist <- 1
        if (x.limDist >= 100)                       x.tickDist <- 10
    
        x.tickDist.5 <- x.tickDist * 5					## for rounding use 5x the distance between the smaller ticks?
        x.tickLim <- c( floor(  min(x)/x.tickDist.5)*x.tickDist.5, 	## round x-axis limits using the distance between ticks
                        ceiling(max(x)/x.tickDist.5)*x.tickDist.5 )
    
    ## plot -- setup plot with title and labels
    plot(0 , xlab=x.lab , ylab = "Hazard Ratio" , main = title , type = "n" , xlim=x.lim , ylim=c(-0.2,y.lim[2]), cex=1.2)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =bg_box_col)
    mtext( subtitle, cex=1.1 )     
    
    # axis additional ticks & reference line
    axis(1, at=seq(x.tickLim[1], x.tickLim[2], x.tickDist), labels=F, tck=-0.009)
    axis(2, at=seq(y.lim[1], y.lim[2], 0.5), labels=F, tck=-0.009)
    axis(2, at=seq(y.lim[1], y.lim[2], 1), labels=F)
	if (y.lim[2] < 5)  axis(2, at=seq(y.lim[1], y.lim[2], 0.25), labels=F, tck=-0.009)
    abline( h = 1 , col = null_line_col , lty = null_line_type , lwd = null_line_wd)

	# adding any custom horizontal or vertical lines?
	if (!is.null(h_lines_pos) & !is.null(h_lines_type) & !is.null(h_lines_col))
	{
		for (hh in 1:length(h_lines_pos))  abline( h = h_lines_pos[hh], col = h_lines_col , lty = h_lines_type , lwd = h_lines_wd)
	}
	if (!is.null(v_lines_pos) & !is.null(v_lines_type) & !is.null(v_lines_col))
	{
		for (vv in 1:length(v_lines_pos))  abline( v = v_lines_pos[vv], col = v_lines_col , lty = v_lines_type , lwd = v_lines_wd)
	}
	
    # add CoxPH results/fitted lines
    lines( sm.spline(x , exp(predicted$fit)) , col = "red" , lwd = 1.5)
    lines( sm.spline(x , exp(predicted$fit + 1.96 * predicted$se)) , col = "orange" , lty = 2 , lwd = 1.5)
    lines( sm.spline(x , exp(predicted$fit - 1.96 * predicted$se)) , col = "orange" , lty = 2 , lwd = 1.5)
    
    # add results legend?
    if (!is.null(legendPos))
    {
      legend(legendPos, inset=0.03,
             paste0("Spline P = ", signif(coef(summary(fit))[1,6] , 2)),
             col="black",
             cex=1, 
             bty="n"
             )
    }
    
    # plot boxplot to show distribution
    if (distBoxPlot)  boxplot(x, horizontal=TRUE, axes=FALSE, pch=19, cex=0.5, add=TRUE, at=-0.2, boxwex=y.lim[2]/20, col=bg_box_col)
    
    # return x and y values?
    if (returnValues) 
    {
      df <- data.frame(x=smsp$x, y=smsp$ysmth, y.u95=smsp.u95$ysmth, y.l95=smsp.l95$ysmth)
      xMinY <- df$x[ df$y == min(df$y) ]
      print( paste0( "X value with minimum HR = ", xMinY ) )      
      return( df )
    }
}

