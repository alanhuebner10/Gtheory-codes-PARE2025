library(gtheory)
text_dat <- read.csv("text_dat_long.csv", header = T)

#The random effects G-theory model presented in Table 1 may be fit with the following code:
formula <- "TextQual ~ (1|id) + (1| scale)+(1| rater) + (1|id: scale) + 
                    (1|id: rater) + (1|rater: scale)"
g_study <- gstudy(data = text_dat, formula)

#View G study results
g_study 
  
#For the D study, specify the sample sizes nr and ns. Here, we use the actual values in data. 
#(They can also be set to a potential sample sizes)
nr <- length(unique(text_dat$rater))
ns <- length(unique(text_dat$scale))
  
#Extract the individual variance components from the model summary
var_p <- g_study$components[g_study$components$source == "id", 2]
var_r <- g_study$components[g_study$components$source == "rater", 2] 
var_s <- g_study$components[g_study$components$source == "scale", 2]
var_pr <- g_study$components[g_study$components$source == "id:rater", 2]
var_ps <- g_study$components[g_study$components$source == "id:scale", 2]
var_rs <- g_study$components[g_study$components$source == "rater:scale", 2]
var_prs <- g_study$components[g_study$components$source == "Residual", 2]

# Then, compute the D study quantities using the variance components and the sample sizes:
var_P <- var_p
var_R <- var_r/nr
var_S <- var_s/ns
var_PR <- var_pr/nr
var_PS <- var_ps/ns
var_RS <- var_rs/(ns*nr)
var_PRS <- var_prs/(ns*nr)

##Finally, compute the reliability coefficients for with s random and fixed:
#Random model (both raters and rating scales fixed)
var_P/(var_P + var_R + var_S + var_PR + var_PS + var_RS + var_PRS)
##Mixed model (rating scales fixed)
(var_P + var_PS)/(var_P + var_R + var_PR + var_PS + var_RS + var_PRS)

  
  