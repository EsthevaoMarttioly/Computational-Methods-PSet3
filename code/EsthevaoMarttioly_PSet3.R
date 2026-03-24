########## Computational Methods - Problem Set 3 - Esthevao Marttioly ##########

rm(list = ls())       ## Be careful with this! It clears the environment

# renv::init()        ## Freeze the package version == just in the first time
# renv::snapshot()    ## Look the package version   == just in the first time
# renv::restore()     ## To restore, Answer "1"     == just in the first time of a new computer


# Set a seed for future replications
set.seed(20260318)


# Import Packages
library(tidyverse)
library(nloptr)
library(ggplot2)
library(latex2exp)
library(stargazer)
library(statmod)


###### Save a theme for the graphs ######
mytheme = theme(legend.position = "bottom",
                plot.title = element_text(size = 12, face = "bold"),
                plot.subtitle = element_text(size = 10),
                panel.background = element_rect(fill = "transparent", colour = "black",
                                                linewidth = 0.5, linetype = "solid"),
                panel.grid.major.y = element_line(colour = "grey", linewidth = 0.5),
                panel.grid.minor.y = element_line(colour = "grey", linewidth = 0.5),
                panel.grid = element_line(colour = "grey98"),
                panel.grid.major.x = element_line(colour = "transparent"),
                panel.grid.minor.x = element_line(colour = "transparent"),
                axis.text = element_text(colour = "black", size = 9),
                strip.background = element_rect(fill = "grey95", colour = "black"),
                strip.text = element_text(colour = "black", size = 9))

colours = c("#85C0F9", "#0F2080", "#A95AA1", "#F5793A", "slateblue2")
# ggsave("output/name.png", width = 5, height = 4)   ## Save the graph into png

#########################################
####### 1. Optimization Methods #########
#########################################

# Defining the function
func1 = function(x) x * sin(5*x)

## Derivatives of the function
Dfunc1 = function(x) sin(5*x) + 5*x * cos(5*x)

D2func1 = function(x) 10*cos(5*x) - 25*x*sin(5*x)


# Understanding visually its behavior
seq(0, 10, length.out = 500) %>%
  data.frame() %>% mutate(y = func1(.)) %>%
  ggplot(aes(x = ., y = y)) + mytheme +
  geom_line(linewidth = 1, colour = colours[2]) +
  geom_hline(yintercept = 0) +
  labs(y = "f(x)", x = "x",
       title = TeX("Function $f(x) = x \\cdot sin(5x), x \\in [0,10]$"))

ggsave("output/fig/p3_q1_function.png", width = 5, height = 4)   ## Save the graph into png


### 1. Grid Search ### - It doesn't have tolerance 

grid_search = function(f, lower, upper, n = 10000){
  
  if(lower >= upper) stop("lower must be < upper")
  
  x_grid = seq(lower, upper, length.out = n)
  f_grid = f(x_grid)
  
  idx = which.min(f_grid) # Index of Min
  
  return(list(x = x_grid[idx],
              fval = f_grid[idx]))
}

#### Run
result1_grid = grid_search(func1, 0, 10)


### 2. Newton Method (from previous PSet) ###

newton_raphson = function(func, deriv_func, x0, tol = 1e-8, max_iter = 1000){
  
  time_start = Sys.time()
  iter = 0
  
  # First Iteration
  x1 = x0 + 10*tol
  converged = F
  
  while (!converged && iter < max_iter) {
    
    fx = func(x0)
    dfx = deriv_func(x0)
    
    if(abs(dfx) < 1e-14) stop("Derivative close to zero")
    
    x1 = x0
    x0 = x0 - fx/dfx      # xn+1 = xn + f(xn) / f'(xn)
    
    iter = iter + 1
    
    # Absolute and Relative Tolerance
    conv_x = abs(x1 - x0) <= max(tol, tol * abs(x1))
    conv_f = abs(deriv_func(x0)) <= max(tol, tol * abs(func(x0)))
    converged = conv_x || conv_f
    
  }
  
  # Check if the solution converged
  if (iter == max_iter) {
    stop(paste("Newton-Raphson method did not converge within", max_iter, "iterations."))
  }
  
  return(list(root = x0, iterations = iter,
              time = as.numeric(Sys.time() - time_start, units = "secs")))
}

#### Run
result1_newt = newton_raphson(Dfunc1, D2func1, 0, tol = 1e-4)


### 3. Controlled Random Search ###

result1_crs = nloptr(x0 = 0, eval_f = func1, lb = 0, ub = 10,
                     opts = list("algorithm" = "NLOPT_GN_CRS2_LM",
                                 "ftol_rel" = 1e-4, "ftol_abs" = 1e-4,
                                 "xtol_rel" = 1e-4, "xtol_abs" = 1e-4,
                                 "nlopt_set_population" = 1000))

### Aggregated Results

cbind(c("Solution", "Objective", "Iterations"),
      c(result1_grid, 1),     ## Grid Search does just 1 iteration
      c(result1_newt$root, func1(result1_newt$root), result1_newt$iterations),
      c(result1_crs$solution, result1_crs$objective, result1_crs$iterations)) %>%
  as.data.frame() %>% `colnames<-`(c("", "Grid Search", "Newton", "CRS")) %>%
  stargazer(summary = F, digits = 4, rownames = F, out = "output/table/q1.tex",
            label = "tab:optmethods", title = "Optimization Methods for $f(x)$")


### Trying another x0 for Newton

newton_raphson(Dfunc1, D2func1, 9, tol = 1e-4)


####################################################
####### 2. Optimization in More Dimensions #########
####################################################

# Defining the function
func2 = function(x1, x2, theta) {
  if (!is.numeric(theta) || length(theta) != 4) {
    stop("theta must be a numeric four-dimension parameter.")
  }
  
  fval = theta[1] * x1 + theta[2] * exp(-x2*x2) + theta[3] * log(1+abs(x1)) + theta[4] * x1^x2
  
  return(fval)
}


# Observable Points
data_points = data.frame(x1 = c(1, 2, -1, 2),
                         x2 = c(1, 4, 2, -2),
                         y  = c(43.614, 563.694, 43.230, 23.130))


# Trace Environment
trace_env = new.env()

trace_env$iter  = 0
trace_env$theta = list()
trace_env$value = numeric()


# Loss Function
loss_func = function(theta, func = func2, data = data_points) {
  if (!is.numeric(theta) || length(theta) != 4) {
    stop("theta must be a numeric four-dimension parameter.")
  }
  
  residuals = func(data[1], data[2], theta) - data_points[3]
  value = sum(residuals^2)
  
  return(value)
}


loss_func_trace = function(theta) {
  value = loss_func(theta)
  
  trace_env$theta[[trace_env$iter + 1]] = theta
  trace_env$value[trace_env$iter + 1] = value
  trace_env$iter = trace_env$iter + 1
  
  return(value)
}


### (b)

loss_func(rep(0, 4))   ## = 322,056.9


### (c) Simulated Annealing (SAMIN) Algorithm

time_start = Sys.time() # Set start time

result2_loss = optim(rep(0, 4), loss_func_trace,
                     method = "SANN", lower = -Inf, upper = Inf,
                     control = list(temp = 20, maxit = 50000))
# control = list(trace = T, REPORT = 1)) would print trace information

time = as.numeric(Sys.time() - time_start, units = "secs") # calculate time
time


### Converting Trace to df

trace_df = cbind.data.frame(1:trace_env$iter,
                            do.call(rbind, trace_env$theta),
                            trace_env$value) %>%
  `colnames<-`(c("iter", paste0(TeX("$theta$"), 1:4), "loss"))


### Trace Graphics

trace_df %>%
  ggplot(aes(x = iter/1e3, y = loss)) +
  geom_line(linewidth = 1, colour = colours[2]) +
  geom_hline(yintercept = 0) + ylim(0, 1000) + mytheme +
  labs(title = "Loss function convergence over iterations",
       x = "Iteration (thousand)", y = TeX("$g(\\theta_i)$"))

ggsave("output/fig/p3_loss_conveg.png", width = 5, height = 4)   ## Save the graph into png


trace_df %>%
  pivot_longer(cols = starts_with("theta")) %>%
  ggplot(aes(x = iter/1e3, y = value, colour = name)) + mytheme +
  geom_line(linewidth = 1) + theme(legend.position = "") + 
  facet_wrap(~name, scales = "free_y") +
  labs(title = TeX("$\\theta_i$'s convergence over iterations"),
       x = "Iteration (thousand)", y = TeX("$(\\theta_i)$")) +
  scale_colour_manual(values = colours[1:4])

ggsave("output/fig/p3_trace_loss.png", width = 5, height = 4)   ## Save the graph into png


trace_df[trace_df$loss == result2_loss$value,]$iter  ## It converged in Iteration 30,195


##########################################
####### 3. Numerical Expectation #########
##########################################

# Gauss Hermite
func3_gh = function(n = 20) {
  gh = gauss.quad(n, kind = "hermite")
  z = gh$nodes
  w = gh$weights
  
  ## Standard Normal
  x = sqrt(2) * z
  y = sqrt(2) * z
  
  ## Double sum
  result3_gh = 0
  
  for (i in 1:n) {
    for (j in 1:n) {
      result3_gh = result3_gh + w[i] * w[j] * max(x[i], y[j])
    }
  }
  
  return(result3_gh / pi) # Normalization
}

# Monte Carlo
func3_mc = function(n = 1e6) {
  X = rnorm(n)
  Y = rnorm(n)
  
  u = pmax(X, Y)
  return(mean(u))
}

# Print results
func3_gh(20)      # 0.5526643
func3_mc(1e6)     # 0.5642427


##########################################
####### 4. Numerical Integration #########
##########################################

# Trapezoid Rule
trapezoid_rule = function(f, a, b, n) {
  
  if (n < 2) stop("n must be >= 2")
  
  x = seq(a, b, length.out = n)
  h = (b - a) / (n - 1)
  
  y = f(x)
  
  result = h * (0.5 * y[1] + sum(y[2:(n-1)]) + 0.5 * y[n])
  
  return(result)
}


# Functions to be approximated
func4_a = function(x) x
func4_b = function(x) x * sin(x)
func4_c = function(x) sqrt(1 - x^2)


# Nodes
n_values = c(3, 5, 10, 15, 20)


# Calculate the Rule
results4 = data.frame(n = n_values,
                      integral_a = NA,
                      integral_b = NA,
                      integral_c = NA)


for (i in seq_along(n_values)) {
  n = n_values[i]
  
  results4$integral_a[i] = trapezoid_rule(func4_a, 0, 1, n)
  results4$integral_b[i] = trapezoid_rule(func4_b, 0, 1, n)
  results4$integral_c[i] = trapezoid_rule(func4_c, 0, 1, n)
}


# Print results
results4 = rbind(c("Analytical", 0.5, sin(1) - cos(1), pi/4), results4) %>%
  `colnames<-`(c("Nodes", "fa", "fb", "fc"))

results4$fb = round(as.numeric(results4$fb), 6)
results4$fc = round(as.numeric(results4$fc), 6)

results4

stargazer(results4, summary = F, digits = 6, rownames = F,
          out = "output/table/q4.tex", label = "tab:trapezoidrule",
          title = "Numerical Integration - Trapezoid Rule")


##############################################
####### 5. Numerical Differentiation #########
##############################################

# Two pointed centered
cd2 = function(f, x, h) (f(x + h) - f(x - h)) / (2*h)

# Four pointed centered
cd4 = function(f, x, h) {
  (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h)
}


# Functions to be approximated
func5_a = function(x) x^2
func5_b = function(x) log(x)
func5_c = function(x) x * sin(x)

x_a = 5
x_b = 10
x_c = 1


# Calculate 
h_values = c(0.001, 0.005, 0.01, 0.05)

results5 = data.frame()

for (h in h_values) {
  
  aux = data.frame(h = paste0("h = ", h),
                   a_cd2 = cd2(func5_a, x_a, h),
                   a_cd4 = cd4(func5_a, x_a, h),
                   b_cd2 = cd2(func5_b, x_b, h),
                   b_cd4 = cd4(func5_b, x_b, h),
                   c_cd2 = cd2(func5_c, x_c, h),
                   c_cd4 = cd4(func5_c, x_c, h))
  
  results5 = rbind(results5, aux)
}


# Print results
results5 = cbind(c("Function", "fa with 2 points", "fa with 4 points",
                   "fb with 2 points", "fb with 4 points",
                   "fc with 2 points", "fc with 4 points"),
                 c("Analytical", rep(c(2 * x_a, 1 / x_b, sin(x_c) + x_c * cos(x_c)), each = 2)),
                 t(results5))

# Clear the output
colnames(results5) = results5[1,]       # Set first row as title
results5 = results5[2:nrow(results5),]  # Delete first row

results5[,2:ncol(results5)] = round(as.numeric(results5[,2:ncol(results5)]), 6)


# Print output
results5

stargazer(results5, summary = F, digits = 4, rownames = F,
          out = "output/table/q5.tex", label = "tab:differentiation",
          title = "Numerical Differentiation")



