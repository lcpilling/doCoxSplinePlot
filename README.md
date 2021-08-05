doCoxSplinePlot
==============

R function to plot the (spline smoothed) exposure in a Cox proportional hazard model against the hazard ratio.

### Description
The function takes as input the results of a Cox proportional hazard model and plots a continuous exposure against the hazard ratio. The continuous exposure must be a spline term for the smoothing function to work. It is up to you to create the sensible CoxPH model.

Includes 95% confidence intervals and boxplot along x-axis to show data distribution.

### Essential inputs
```
# x        vector of x-values (the risk factor)
# fit      CoxPH object -- fitted model of x against mortality
# x.lab    label for x-axis, e.g. "Serum albumin (g/dL)"
# title    title for the plot
```

### Optional inputs for customisation of plot
```
# subtitle        additional subtitle
# returnValues    return a data.frame of the x and y values (default is FALSE)
# legendPos       default is NO results legend. Specify "topleft" or similar to show the results on the plot
# distBoxPlot     show bowplot of distribution? Default is FALSE
# yMax            manually set the y-axis maximum? (minimum always 0)
# xLims           manually set the x-axis limits? (needs to be a vector of length 2 - min or max can be NULL)
# bg_box_col      specific background colour of plot (default = grey98)
# null_line_type  horizontal reference line type (default = 2)
# null_line_col   horizontal reference line colour (default = grey60)
# null_line_wd    horizontal reference line width (default = 1.5)
# h_lines_pos     vector of y axis positions for grid lines if desired (default = NULL)
# h_lines_type    type of lines for y axis (default = 2)
# h_lines_col     colour of lines for y axis (default = grey60)
# h_lines_wd      width of lines for y axis (default = 1)
# v_lines_pos     vector of x axis positions for grid lines if desired (default = NULL)
# v_lines_type    type of lines for x axis (default = 2)
# v_lines_col     colour of lines for x axis (default = grey60)
# v_lines_wd      width of lines for x axis (default = 1)
# ...             other plot options (e.g. cex)
```

### Depends upon packages
- survival
- pspline

### Example
```
## get data, exclude missings
data <- na.omit( data.frame( dead,       # numeric vector: binary [0=alive, 1=dead]
                             age_death,  # numeric vector: age at death for each participant
                             albumin     # numeric vector: the risk factor you are assessing
                             age,        # numeric vector: age at measurement of risk factor (cofactor)
                             sex,        # numeric vector: sex of participant (cofactor)
                             smokes      # numeric vector: smoking status of participant (cofactor)
                           ) )

## create the 'survival object' -- load survival package
require(survival)
surv.death <- Surv(data$age_death , data$dead)

## do Cox model -- can include cofactors if desired
##     primary independent variable (exposure of interest) must be pspline()
pham.fit <- coxph( surv.death ~ pspline(data$albumin, df=4) 
                                + data$age + data$sex + as.factor(data$smokes) )

## use doCoxSplinePlot() to plot the smoothed curve (including CI's) for the Cox model
doCoxSplinePlot(x        = data$albumin, 
                fit      = pham.fit, 
                x.lab    = "Albumin (g/dL)", 
                title    = "Circulating albumin and mortality risk",
                subtitle = "CoxPH model adjusted for age, sex and smokes, with 95% CIs")

```
![](https://github.com/lukepilling/doCoxSplinePlot/blob/master/I0wQG.png)
