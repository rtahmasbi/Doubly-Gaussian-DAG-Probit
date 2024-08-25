################################
# Rscript simulation.R 30 1000 500 .1 3000 1000 .1
################################
print('==============================================================================')
print('> start ...')

library(pcalg)
library(mvtnorm)
library(miscTools)
library(resample)
library(DensParcorr)


args  = commandArgs(trailingOnly = TRUE)
print(args)
q = as.numeric(args[1])
n1 = as.numeric(args[2])
n2 = as.numeric(args[3])
prob_edge = as.numeric(args[4])
TT = as.numeric(args[5])
burn = as.numeric(args[6])
cor_a1_a2 = as.numeric(args[7])

if (FALSE) {
	q = 30   # number of nodes
	n1 = 1000 # sample size
	n2 = 1000 # sample size
	prob_edge = .1 # 6/(2*q-2)
	TT   = 10000
	burn = 4000
	cor_a1_a2 = .1
}


rand_int = sample.int(1e6, 1, replace = TRUE)

a = q
g1 = 1/n1
g2 = 1/n2


# Generate true DAG and parameters

#set.seed(1)

make_A = function(q, prob) {
	true_dag = randomDAG(q, prob = prob)
	A = t(as(true_dag, "matrix"))
	A[A != 0] = 1
	return(A)
}


lower_elements = function(M) {
	ret = M[lower.tri(M)][order(-row(M)[lower.tri(row(M))])]
	return (ret)

}

diag_elements = function(M) {
	ret = M[row(M)==col(M)]
	return (ret)

}


pa = function(set, object) {
  pa = which(object[,set] != 0)
  return(pa)
}

print('> making A1 ...')
A1 = make_A(q, prob= prob_edge)
print('> making A1 done')

print('> making A2 ...')
# make correlated A's
repeat {
	A2 = make_A(q, prob=prob_edge)
	cr = cor(lower_elements(A1),lower_elements(A2))
	if (cr>cor_a1_a2) {
		break
	}
}
print('> making A2 done')


temp_L = matrix(runif(q*q, 1, 2), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE)
L1 = A1* temp_L; diag(L1) = 1
L2 = A2* temp_L; diag(L2) = 1
D = diag(rep(1, q))
#D = diag(runif(q, 1, 2))

Sigma1 = solve(t(L1))%*%D%*%solve(L1)
Sigma2 = solve(t(L2))%*%D%*%solve(L2)
Omega1 = solve(Sigma1)
Omega2 = solve(Sigma2)

true_Partial1 = prec2part(Omega1)
true_Partial2 = prec2part(Omega2)



# normal conditional distribution
Sigma_bar = function(j, A, Sigma)
{
	pa_j = pa(j, A)
	aa = matrix(Sigma[pa_j,j], nrow=1)
	return(Sigma[j,j] - aa %*% solve(Sigma[pa_j,pa_j]) %*% t(aa))
}



mu1    = c(rep(0, q))
mu2    = c(rep(0, q))

# Generate the data


Y1 = rmvnorm(n1, mu1, Sigma1)
Y2 = rmvnorm(n2, mu2, Sigma2)

m1 = colMeans(Y1)
m2 = colMeans(Y2)
#s = apply(X = Y, FUN = sd, MARGIN = 2)

Y1 = t(t(Y1) - m1)
Y2 = t(t(Y2) - m2)

#head(Y1)


#################
## MCMC scheme ##
#################

## Finally, create the 0-1 response variable by thresholding of Y[,1]

X1 = Y1[,-1]
X2 = Y2[,-1]
y1 = Y1[,1]
y2 = Y2[,1]
y1[y1 > 0] = 1
y1[y1 < 0] = 0
y2[y2 > 0] = 1
y2[y2 < 0] = 0


## Data consists of y (binary response) and X (n,p) matrix with covariates


## Fix number of MCMC iterations and burn-in period


source("mcmc_dag_probit.R")

## Posterior inference on DAGs and parameters (Sigma)

#t_0 = proc.time()
#out = mcmc_dag_probit(y1 = y1, y2 = y2, X1 = X1, X2 = X2, TT = TT, burn = burn, causal = FALSE)
#t_1 = proc.time() - t_0
#t_1

## Posterior inference on DAGs and parameters (Sigma) + causal effect estimation

gamma_0 = runif(1)*2-1

t_0 = proc.time()
out_causal = mcmc_dag_probit(y1 = y1, y2 = y2, X1 = X1, X2 = X2, TT = TT, burn = burn, gamma_0= gamma_0)
t_1 = proc.time() - t_0
print(t_1)

#################################
## Compute posterior summaries ##
#################################

## Posterior probabilities of edge inclusion

probs1 = round(matrix(rowMeans(out_causal$A1_chain[,(burn + 1):TT]), q, q), 4)
probs2 = round(matrix(rowMeans(out_causal$A2_chain[,(burn + 1):TT]), q, q), 4)


## Median Probability DAG Model

A1_hat = round(probs1)
A2_hat = round(probs2)

print('> A1_hat')
print(mean(A1_hat[lower.tri(A1_hat, diag = FALSE)] == A1[lower.tri(A1, diag = FALSE)]))
print('> A2_hat')
print(mean(A2_hat[lower.tri(A2_hat, diag = FALSE)] == A2[lower.tri(A2, diag = FALSE)]))


## BMA estimate of the covariance matrix

Sigma1_hat = matrix(rowMeans(out_causal$Sigma1_chain[,(burn + 1):TT]), q, q)
Sigma2_hat = matrix(rowMeans(out_causal$Sigma2_chain[,(burn + 1):TT]), q, q)


Omega1_hat = solve(Sigma1_hat)
Omega2_hat = solve(Sigma2_hat)


print('> Sigma1_hat')
print(mean(Sigma1_hat-Sigma1))
print(mean(abs(Sigma1_hat-Sigma1)))



print('> D_chain - rowMeans')
print(round(diag_elements(matrix(rowMeans(out_causal$D_chain[,(burn + 1):TT]), q, q)), 2))
print('> D_chain - rowMedians')
print(round(diag_elements(matrix(rowMedians(out_causal$D_chain[,(burn + 1):TT]), q, q)), 2))




## BMA estimate of causal effects

#print(round(out_causal$Causal1_hat, 2))

print('> Causal1_hat')
print(round(colMeans(out_causal$Causal1_hat), 2))
print('> Causal2_hat')
print(round(colMeans(out_causal$Causal2_hat), 2))


print('> Gamma_chain')
print(mean(out_causal$Gamma_chain[(burn + 1):TT]))



L_check = function(L, L_hat, A) {
	x_true = L[A==1]
	x_pred = L_hat[A==1]
	y_true = L[lower.tri(A, diag = FALSE) & A1==0]
	y_pred = L_hat[lower.tri(A, diag = FALSE) & A1==0]
	a1 = mean(x_true-x_pred)
	a2 = var(x_true-x_pred)
	a3 = mean(abs(x_true-x_pred))
	b1 = mean(y_true-y_pred)
	b2 = var(y_true-y_pred)
	b3 = mean(abs(y_true-y_pred))
	return(c(a1,a2,a3,b1,b2,b3))
}

c_general = data.frame(t(c(rand_int=rand_int , q=q, n1=n1, n2=n2, prob_edge=prob_edge, TT=TT, burn=burn, cor_a1_a2=cor_a1_a2)))
#colnames(c_general) = c("rand_int" , "q", "n1", "n2", "prob_edge", "TT", "burn")

#print(round(matrix(rowMeans(out_causal$L1_chain[,(burn + 1):TT]), q, q), 2))

L1_hat = matrix(rowMeans(out_causal$L1_chain[,(burn + 1):TT]), q, q)
L2_hat = matrix(rowMeans(out_causal$L2_chain[,(burn + 1):TT]), q, q)
#mean(abs(L1_hat[lower.tri(L1_hat, diag = FALSE)] == A1[lower.tri(L1, diag = FALSE)]))
res_L1 = L_check(L1, L1_hat, A1)
res_L2 = L_check(L2, L2_hat, A2)
RES_L = data.frame(t(c(res_L1, res_L2)))


tbl_conf1 = table(Predictions = A1_hat[lower.tri(A1_hat, diag = FALSE)], TrueLabels = A1[lower.tri(A1, diag = FALSE)])
tbl_conf2 = table(Predictions = A2_hat[lower.tri(A2_hat, diag = FALSE)], TrueLabels = A2[lower.tri(A2, diag = FALSE)])

RES_time = data.frame(t(c(t_1)))

# tbl_conf1[2,2]: TP  X1
# tbl_conf1[1,1]: TN  X2
# tbl_conf1[2,1]: FP  X3
# tbl_conf1[1,2]: FN  X4
RES_conf = data.frame(t(c(TN1=tbl_conf1[1,1], FN1=tbl_conf1[1,2], FP1=tbl_conf1[2,1], TP1=tbl_conf1[2,2],
						  TN2=tbl_conf2[1,1], FN2=tbl_conf2[1,2], FP2=tbl_conf2[2,1], TP2=tbl_conf2[2,2])))

# Sigma
a1 = mean(Sigma1_hat-Sigma1)
a2 = mean(abs(Sigma1_hat-Sigma1))
a3 = var(as.vector(Sigma1_hat-Sigma1))
b1 = mean(Sigma2_hat-Sigma2)
b2 = mean(abs(Sigma2_hat-Sigma2))
b3 = var(as.vector(Sigma2_hat-Sigma2))
RES_sigma = data.frame(t(c(mean_diff1=a1, mean_abs_diff1=a2, var_diff1=a3, mean_diff2=b1, mean_abs_diff2=b2, var_diff2=b3)))


# Omega
a1 = mean(Sigma1_hat-Omega1)
a2 = mean(abs(Sigma1_hat-Omega1))
a3 = var(as.vector(Sigma1_hat-Omega1))
b1 = mean(Sigma2_hat-Omega2)
b2 = mean(abs(Sigma2_hat-Omega2))
b3 = var(as.vector(Sigma2_hat-Omega2))
RES_Omega = data.frame(t(c(mean_diff1=a1, mean_abs_diff1=a2, var_diff1=a3, mean_diff2=b1, mean_abs_diff2=b2, var_diff2=b3)))


#a1 = colMeans(out_causal$Causal1_hat)
#b1 = colMeans(out_causal$Causal2_hat)
#RES_Causal = data.frame(Causal1_hat=a1, Causal2_hat=b1)
causal1_true = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma1, A = A1, X_mat = X1, g = 0))
a1 = colVars(causal1_true)
a2 = colVars(out_causal$Causal1_hat)
causal2_true = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma2, A = A2, X_mat = X2, g = 0))
b1 = colVars(causal2_true)
b2 = colVars(out_causal$Causal2_hat)
RES_Causal = data.frame(A1_c1 = A1[,1], Causal1_true=a1, Causal1_hat=a2, A2_c1 = A2[,1], Causal2_true=b1, Causal2_hat=b2)




a1 = mean(out_causal$Gamma_chain[(burn + 1):TT])
a2 = median(out_causal$Gamma_chain[(burn + 1):TT])
a3 = var(out_causal$Gamma_chain[(burn + 1):TT])
RES_gamma = data.frame(t(c(gamma_0 , a1, a2, a3)))


cr = cor(lower_elements(A1),lower_elements(A2))
RES_A1_A2_corr = data.frame(cr=t(cr))


P1_post = matrix(rowMedians(out_causal$Partial1_chain[,(burn + 1):TT]), q, q)
P2_post = matrix(rowMedians(out_causal$Partial2_chain[,(burn + 1):TT]), q, q)
a1 = mean(lower_elements(true_Partial1) - lower_elements(P1_post))
a2 = mean(abs(lower_elements(true_Partial1) - lower_elements(P1_post)))
a3 = var(lower_elements(true_Partial1) - lower_elements(P1_post))
b1 = mean(lower_elements(true_Partial2) - lower_elements(P2_post))
b2 = mean(abs(lower_elements(true_Partial2) - lower_elements(P2_post)))
b3 = var(lower_elements(true_Partial2) - lower_elements(P2_post))
RES_Partial = data.frame(t(c('mean1'=a1, 'mean_abs1'=a2, 'var1'=a3, 'mean2'=b1, 'mean_abs2'=b2, 'var2'=b3)))

RES_ROC = data.frame('probs1'=lower_elements(probs1), 'A1'=lower_elements(A1), 'probs2'=lower_elements(probs2), 'A2'=lower_elements(A2) )


fname_out = "RES_time.txt"
write.table(data.frame(c(c_general, RES_time)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_L.txt"
write.table(data.frame(c(c_general, RES_L)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_conf.txt"
write.table(data.frame(c(c_general, RES_conf)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_sigma.txt"
write.table(data.frame(c(c_general, RES_sigma)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_Omega.txt"
write.table(data.frame(c(c_general, RES_Omega)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_Causal.txt"
write.table(data.frame(c(c_general, RES_Causal)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_gamma.txt"
write.table(data.frame(c(c_general, RES_gamma)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_A1_A2_corr.txt"
write.table(data.frame(c(c_general, RES_A1_A2_corr)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_Partial.txt"
write.table(data.frame(c(c_general, RES_Partial)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )

fname_out = "RES_ROC.txt"
write.table(data.frame(c(c_general, RES_ROC)), file = fname_out, append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(fname_out) )


print('> end.')


