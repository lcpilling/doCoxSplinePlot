doCoxSplinePlot
==============

R function to plot the (spline smoothed) exposure in a Cox proportional hazard model against the hazard ratio.

### Description
The function takes as input the results of a Cox proportional hazard model and plots a continuous exposure against the hazard ratio. The continuous exposure must be a spline term for the smoothing function to work. It is up to you to create the sensible CoxPH model.

Includes 95% confidence intervals and boxplot along x-axis to show data distribution.

### Inputs
* x        :: {REQUIRED} vector of x-values (the primary exposure in the Cox model -- used with pspline() function)
* fit		   :: {REQUIRED} CoxPH object -- fitted model
* x.lab		 :: {REQUIRED} label for x-axis, e.g. "Serum albumin (g/dL)" - whatever your exposure is
* title    :: {REQUIRED} title for the plot
* subtitle :: {optional} additional subtitle

### Example
```
## get data, exclude missings
data <- na.omit( data.frame( dead,       # numeric vector: binary [0=alive, 1=dead]
                             age_death,  # numeric vector: age at death for each participant
                             albumin     # numeric vector: the risk factor you are assessing
                             age,        # numeric vector: age at measurement of risk factor (optional cofactor)
                             sex,        # numeric vector: sex of participant (optional cofactor)
                             smokes      # numeric vector: smoking status of participant (optional cofactor)
                           ) )

## create the 'survival object'
surv.death <- Surv(data$age_death , data$dead)

## do Cox model -- can include cofactors if desired -- primary independent variable (exposure of interest) must be pspline()
pham.fit <- coxph( surv.death ~ pspline(data$albumin, df=4) + data$age + data$sex + as.factor(data$smokes) )

## use doCoxSplinePlot() to plot the smoothed curve (including CI's) for the Cox model
doCoxSplinePlot(x=data$albumin, 
                fit=pham.fit, 
                x.lab="Albumin (g/dL)", 
                title="Circulating albumin and mortality risk",
                subtitle="CoxPH model adjusted for age, sex and smokes, with 95% CIs")

```
![](http://s22.postimg.org/vr887q00x/Albumin_mortality_risk.png)
