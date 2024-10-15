library(Rcapture)

# Define a 2x2 contingency table with detected and undetected individuals
# Observed counts for each combination of method detection:
A <- 25404 # detected by Swab
B <- 24440 # detected by N-antibody
C <- 18128 # detected by both
U <- A+B-C # detected by any

data <- data.frame(
  Swab = c(1, 1, 0),         # Detected by swab (1: Yes, 0: No)
  Antibody = c(1, 0, 1),     # Detected by N-antibody (1: Yes, 0: No)
  Freq = c(C, A-C, B-C)       # Frequency: Detected by both, only swab, only N-antibody
)

# Fit a log-linear model to account for potential dependency between methods
cr_model <- profileCI(data, dfreq = TRUE, m = 'Mt')

# "Mt" indicates time or method-dependent model
# using this as expect some heterogeneity in the detection probability,
# e.g. PCR might be better in some than antibody and vice versa. 

# View model output
cr_model

# true infections, point estimate: 
cr_model$results[1]

# so missed infections: 
1-(U / cr_model$results[1])
#CI lower
1-(U / cr_model$results[2])
#CI upper
1-(U / cr_model$results[3])

# missed by swab:
1- (A / cr_model$results[1])
#CI
1- (A / cr_model$results[2])
1- (A / cr_model$results[3])

# missed by antibodies:
1- (B / cr_model$results[1])
#CI
1- (B / cr_model$results[2])
1- (B / cr_model$results[3])
