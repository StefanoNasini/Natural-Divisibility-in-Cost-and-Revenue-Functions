

rm(list=ls(all=TRUE))
gc(reset = TRUE)


library("dplyr")
library(Rmpfr)
library(lpSolve)
library(Matrix)
library(R.utils)

# start /B /HIGH /WAIT Rscript Run_estimation_Banking.R > Run_estimation_Banking.out

###############################################
###############################################


#---------------------------------------------
# Functions
#---------------------------------------------
# Generalized data structure for Cost min and Rev max
#---------------------------------------------

#' Draw a filled polygon as a step function
#'
#' @param x        x coordinates; must be one longer than `y1` and `y2`
#' @param y1       lower y coordinates
#' @param y2       upper y coordinates
#' @param border   Color to draw the border; see `polygon()`.  Default: no borders.
#' @param ...      Additional arguments to be passed to `polygon()`
#'
#' @examples
#' plot(1:10, 1:10, type='n', ylim=c(0,12))
#' polygon.step(1:10, 0:8, 2:10, col='gray')
#' lines(1:10, 1:10, type='s')
#'
#' @export

polygon.step <- function(x, y1, y2, border=FALSE, ...) {
  nx <- length(x)
  ny <- length(y1)
  if (length(y2)!=ny) stop("y1 and y2 must be the same length")
  if (nx != (ny+1)) stop("x must be one longer than y")
  xx <- c(x[1], rep(x[-c(1,nx)], rep(2,nx-2)), x[nx])
  xxx <- c(xx, rev(xx))
  yy1 <- rep(y1, rep(2,ny))
  yy2 <- rep(y2, rep(2,ny))
  yyy <- c(yy1, rev(yy2))
  polygon(xxx, yyy, border=border, ...)

}


A_mat = function(s){

	if(s >= 1){
		out = t(YY)
	}else{
		out = -t(XX)
	}

	return(out)

}

q_mat = function(s){

	if(s >= 1){
		out = m
	}else{
		out = n
	}

	return(out)

}

p_mat = function(s,kkk){

	if(s >= 1){
		out = as.numeric(P_I[kkk,])
	}else{
		out = -as.numeric(P_O[kkk,])
	}

	return(out)

}

D_mat = function(s){

	if(s >= 1){
		out = diag(n)
	}else{
		out = -diag(m)
	}

	return(out)

}


b_mat = function(s, kkk){

	if(s >= 1){
		out = as.numeric(YY[kkk,])
	}else{
		out = -as.numeric(XX[kkk,])
	}

	return(out)

}



LP_mat = function(s, MM, dd){

	out = rbind( 
		cbind(A_mat(s)*(dd/MM), matrix(0, q_mat(s), q_mat(1-s))),
		cbind(A_mat(1-s)*(dd/MM), D_mat(s)),
		c(rep(1, k), rep(0, q_mat(1-s))),
		cbind(diag(rep(1, k)), matrix(0, k, q_mat(1-s)))
	)

	return(out)

}

LP_mat_eq = function(s, MM, dd, slack){

	if(slack){

		out = rbind( 
			cbind(A_mat(s)*(dd/MM), -diag(q_mat(s))),
			c(rep(1, k), rep(0, q_mat(s)))
		)

	}else{

		out = rbind( 
			cbind(A_mat(s)*(dd/MM)),
			c(rep(1, k))
		)

	}

	return(out)

}


#---------------------------------------------
# Compute utilities
#---------------------------------------------


Relaxed_Efficiency = function(ddd, MMM, kk, ss){

	LLP_matrix = LP_mat(ss, MMM, ddd) 

	cc_k = c(rep(0,k), p_mat(ss,kk))
	ccc_k = c(cc_k, rep(0,k + q_mat(1-ss) ))

	LLLP_matrix = bdiag(LLP_matrix, LLP_matrix)

	LLLP_matrix = as.matrix(rbind(LLLP_matrix, LL_second_block, LL_second_block  ))
	
	#-----------------------------------------------


	if(MMM <= 20){

		bb_k = c(b_mat(ss,kk), rep(0, q_mat(1-ss)), rep(MMM, k + 1))	
		bbb_k = c(bb_k, bb_k, aa*rep(1,k), -aa*rep(1,k) )		

		out_EDM = lp(direction = "min", objective.in = ccc_k, const.mat = LLLP_matrix, const.dir = DIR_M_relax, const.rhs = bbb_k, int.vec = INTEGER_var)

		if(out_EDM$status < 2){

			out1 = sum(
				out_EDM$solution[INTEGER_var]* (
					c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))
				))

		}else{
			out1  = 1.0e+20
		}

	}else{

		LLP_matrix_M = LP_mat_eq(ss, 1, ddd, slack = FALSE) 

		bb_k = c(b_mat(ss,kk), 1)

		cc_k = c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))

		out_DEA = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix_M, const.dir = DIR, const.rhs = bb_k)

		if(out_DEA$status < 2){
			out1  = as.numeric(out_DEA$objval) 
		}else{
			out1  = 1.0e+20
		}

	}

	return(out1/MMM)

}



Efficiency = function(ddd, MMM, kk, ss){

	LLP_matrix = LP_mat(ss, MMM, ddd) 
	bb_k = c(b_mat(ss,kk), rep(0, q_mat(1-ss)), rep(MMM, k + 1))
	cc_k = c(rep(0,k), p_mat(ss,kk))

	#-----------------------------------------------

	if(MMM <= 20){		

		out_EDM = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix, const.dir = DIR_M, const.rhs = bb_k, int.vec = 1:k)
		
		if(out_EDM$status < 2){
			out1  = as.numeric(out_EDM$objval) 
		}else{
			out1  = 1.0e+20
		}

	}else{		

		LLP_matrix_M = LP_mat_eq(ss, 1, ddd, slack = FALSE) 

		bb_k = c(b_mat(ss,kk), 1)

		cc_k = c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))

		out_DEA = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix_M, const.dir = DIR, const.rhs = bb_k)

		#print(length(bb_k))
		#print(nrow(LLP_matrix_M1))
		#print(length(DIR))
		#print(length(cc_k))
		#print(ncol(LLP_matrix_M1))

		if(out_DEA$status < 2){
			out1  = as.numeric(out_DEA$objval) 
		}else{
			out1  = 1.0e+20
		}

	}

	return(out1)

}


utilities = function(MM, ss, case, PPPi, d_val){

	out1 = 0
	out1_Rel = 0

	print("--------------------------------")
	print(MM)
	print("--------------------------------")
	
	for (kkk in 1:k){

		gc(reset = TRUE)

		if (case == 0){

			if(ss == 0){

				pb = as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				tilde_a_k = abs(PPPi[[kkk]])

				#Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				#Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				tryCatch({
					Pi <- apply_with_timeout(d_val, 0, MM, kkk, ss) #  Efficiency
					# print(c("Good", Pi))

				}, TimeoutException = function(ex) {
					# Assign zero to Pi if the function execution is interrupted due to timeout
					Pi <- 0.5*tilde_a_k + 0.5*pb
					# print(c("Bad", Pi))
				})
				tryCatch({
					Pi_rel <- apply_with_timeout(d_val, 1, MM, kkk, ss)#Relaxed_Efficiency
					# print(c("Good", Pi_rel))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi_rel if the function execution is interrupted due to timeout
					Pi_rel <- 0.5*tilde_a_k + 0.5*pb
					# print(c("bad", Pi_rel))
				})

				pb_Pi = max(0, (pb - min(Pi)))
				pb_Pi_rel = max(0, (pb - min(Pi_rel)))

				OUT = max(1 - pmax(0, 
						pmin( 1, pb_Pi/tilde_a_k  )
					))
				OUT_rel = max(1 - pmax(0, 
						pmin(1, pb_Pi_rel/tilde_a_k, pb_Pi/tilde_a_k )
					))
			
				#print(c(OUT, OUT_rel ))

			}else{

				tilde_a_k = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				pb = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				#Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				#Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				tryCatch({
					Pi <- apply_with_timeout(d_val, 0, MM, kkk, ss) #  Efficiency
					# print(c("Good Pi: ", Pi))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi if the function execution is interrupted due to timeout
					Pi <- 0.5*tilde_a_k + 0.5*pb
					# print(c("Bad Pi:", Pi))
				})
				tryCatch({
					Pi_rel <- apply_with_timeout(d_val, 1, MM, kkk, ss)#Relaxed_Efficiency
					# print(c("Good Pi_rel:", Pi_rel))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi_rel if the function execution is interrupted due to timeout
					Pi_rel <- 0.5*tilde_a_k + 0.5*pb
					# print(c("Bad Pi_rel", Pi_rel))
				})
				
				pb_Pi = max(0, (pb - min(Pi)))
				pb_Pi_rel = max(0, (pb - min(Pi_rel)))

				# print(c("Final", pb_Pi))
				# print(c("Final", pb_Pi_rel))

				OUT = max(1 - pmax(0, 
						pmin( 1, pb_Pi/abs(tilde_a_k)  )
					))
				OUT_rel = max(1 - pmax(0, 
						pmin(1, pb_Pi_rel/abs(tilde_a_k), pb_Pi/abs(tilde_a_k) )
					))
			
				#print(c(OUT, OUT_rel))

			}

			out1  = out1 + OUT

			out1_Rel  = out1_Rel + OUT_rel


		}else{

			if(ss == 0){

				tilde_a_k = abs(PPPi[[kkk]])

				pb = as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				#Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				#Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				tryCatch({
					Pi <- apply_with_timeout(d_val, 0, MM, kkk, ss) #  Efficiency
					# print(c("Good", Pi))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi if the function execution is interrupted due to timeout
					Pi <- 0.5*tilde_a_k + 0.5*pb
					# print(c("bad", Pi))
				})
				tryCatch({
					Pi_rel <- apply_with_timeout(d_val, 1, MM, kkk, ss)#Relaxed_Efficiency
					# print(c("Good", Pi_rel))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi_rel if the function execution is interrupted due to timeout
					Pi_rel <- 0.5*tilde_a_k + 0.5*pb
					# print(c("bad", Pi_rel))
				})
				
				pb_Pi = pmax(0, (pb - Pi))
				pb_Pi_rel = pmax(0, (pb - Pi_rel))

				OUT = max(1 - pmin(0.99999999, pb_Pi/abs(tilde_a_k)))
				OUT_rel = max(1 - pmin(0.99999999, pb_Pi_rel/abs(tilde_a_k)))
			
				#print(c(OUT, OUT_rel ))

			}else{

				tilde_a_k = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))
				pb = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				#Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				#Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				tryCatch({
					Pi <- apply_with_timeout(d_val, 0, MM, kkk, ss) #  Efficiency
					# print(c("Good", Pi))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi if the function execution is interrupted due to timeout
					Pi <- 0.5*tilde_a_k + 0.5*pb
					# print(c("bad", Pi))
				})
				tryCatch({
					Pi_rel <- apply_with_timeout(d_val, 1, MM, kkk, ss)#Relaxed_Efficiency
					# print(c("Good", Pi_rel))
				}, TimeoutException = function(ex) {
					# Assign zero to Pi_rel if the function execution is interrupted due to timeout
					Pi_rel <- 0.5*tilde_a_k + 0.5*pb
					# print(c("bad", Pi_rel))
				})

				pb_Pi = pmax(0, (pb - Pi))
				pb_Pi_rel = pmax(0, (pb - Pi_rel))

				OUT = max(1 - pmin(0.99999999, pb_Pi/abs(tilde_a_k)))
				OUT_rel = max(1 - pmin(0.99999999, pb_Pi_rel/abs(tilde_a_k)))
			
				#print(c(OUT, OUT_rel))
			}

			out1  = out1 + OUT

			out1_Rel  = out1_Rel + OUT_rel

		}
		
				print(c(kkk, out1/kkk))
				
	}

	out2 = 1 - (1/MM)
	print(MM)

	gc()

	return(c(out1/k, out1_Rel/k, out2))

}



time_limit <- 1000

# Define a function that wraps the gmp::apply function with a time limit
apply_with_timeout <- function(d_val, fun, MM, kkk, ss){
	
	if(fun == 0)
	{
		out = withTimeout(gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss), timeout = 10000, onTimeout = "error")
	}
	if(fun == 1)
	{
		out = withTimeout(gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss), timeout = 10000, onTimeout = "error")
	}
	
	return(out)

}




###################################################################
###################################################################
#
# Banking data
#
###################################################################
###################################################################

setwd("Banking_Data")

data = read.table("Aly,Grabowski,Pasurka&Rangan-1990REStat-data.txt", header = FALSE)

Names = c('ID', 'BankID', 'X_1', 'X_2', 'X_3', 'P_1', 'P_2', 'P_3', 'Y_1', 'Y_2', 'Y_3', 'Y_4', 'Y_5', 'branch-dummy')

length(Names)

names(data) = Names

XX = data[,3:5]
YY = data[,9:13]
P_I = data[,6:8] 

m = ncol(YY)
n = ncol(XX)
k = nrow(XX)

kk = 10

k = nrow(XX)
m = ncol(YY)
n = ncol(XX)	

M_max = 20
d_lev = 2
sss = 1

aa = 0.3

DIR_M = c(rep(">=", n + m), "=", rep("<=", k))
DIR = c(rep(">=", q_mat(sss)), "<=")
DIR_M_relax = c(DIR_M, DIR_M, rep("<=", k), rep(">=", k) )
INTEGER_var = (k + q_mat(1-sss) + 1):(k + q_mat(1-sss) + k)

LL_second_block0 = cbind(diag(k), matrix(0,k,q_mat(1-sss)) )
LL_second_block1 = cbind(-diag(k), matrix(0,k,q_mat(1-sss)) )

LL_second_block = cbind(LL_second_block0, LL_second_block1)

d_val = as.matrix(c(1, seq(0.90, 1.1, length = d_lev)))
M_val = as.matrix(1:M_max)

setwd("Results/Banking")

PPPi = list()
PPPi_min = NULL
Pb = NULL
for(uuu in 1:k){
	
	PPPi[[uuu]] = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = 1000, kk = uuu, ss = sss)
	PPPi_min = c(PPPi_min, min(PPPi[[uuu]]))
	Pb = c(Pb, abs(as.numeric(p_mat(sss,uuu)%*%b_mat(1-sss,uuu))))
}


U = matrix(0, M_max, 3)

for(uuu in 1:M_max){
	U[uuu,] = utilities(MM = uuu, ss = 1, case = 1, PPPi = PPPi, d_val = d_val)
	gc(reset = TRUE)
	Sys.sleep(5) 
	gc(reset = TRUE) 
	save(U, file = "Full_output_cost_min_case_1_Banking_5.RData")
}

#----------------------------------------------------

setwd("Results/Banking")

load("Full_output_cost_min_case_0_Banking.RData")

Relaxed_tt = pmin(
	(U[3,]-min(U[3,]))/(max(U[3,]) - min(U[3,])) ,
	(U[2,]-min(U[1,]))/(max(U[1,]) - min(U[1,]))
)

tt = pmin(
	(U[3,]-min(U[3,]))/(max(U[3,]) - min(U[3,])) ,
	(U[1,]-min(U[1,]))/(max(U[1,]) - min(U[1,]))
)

M_Relaxed = which.max(Relaxed_tt)
M_opt = which.max(tt)

t_Relaxed = max(Relaxed_tt)
t_opt = max(tt)

Table_cost = c(M_opt, t_opt, M_Relaxed, t_Relaxed,
							U[1,M_opt], U[2,M_Relaxed], U[3,M_opt], U[3,M_Relaxed],
							max(U[3,]), max(U[1,]),
							min(U[3,]), min(U[1,])
							)
UU = list("U" = U, "Table" = Table_cost)

write.table(Table_cost, file = "Table_cost_case_0_Banking.log")

U_year = t(U)

filenamePLOT = paste("Plot_s_1_case_0_Banking",  ".pdf", sep = "")

pdf(filenamePLOT)

plot(U_year[3,], U_year[1,], type = 's', , col='gray', lwd = 1, xlab = 'Divisibility player', ylab = 'Efficiency player', ylim = c(min(U_year[1,]), max(U_year[2,])) )

polygon.step(U_year[3,], U_year[1,1:19], rep(min(U_year[1,]),19), col='gray')

lines(U_year[3,], U_year[2,], type = 's', col = "black", lty = 4, lwd = 4, xlab = 'Divisibility player', ylab = 'Efficiency player')

lines(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 2 , col = 2)

points(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 6 , col = 2)

dev.off()



#----------------------------------------------------

load("Full_output_cost_min_case_1_Banking.RData")

Relaxed_tt = pmin(
	(U[3,]-min(U[3,]))/(max(U[3,]) - min(U[3,])) ,
	(U[2,]-min(U[1,]))/(max(U[1,]) - min(U[1,]))
)

tt = pmin(
	(U[3,]-min(U[3,]))/(max(U[3,]) - min(U[3,])) ,
	(U[1,]-min(U[1,]))/(max(U[1,]) - min(U[1,]))
)

M_Relaxed = which.max(Relaxed_tt)
M_opt = which.max(tt)

t_Relaxed = max(Relaxed_tt)
t_opt = max(tt)

Table_cost = c(M_opt, t_opt, M_Relaxed, t_Relaxed,
							U[1,M_opt], U[2,M_Relaxed], U[3,M_opt], U[3,M_Relaxed],
							max(U[3,]), max(U[1,]),
							min(U[3,]), min(U[1,])
							)
UU = list("U" = U, "Table" = Table_cost)

write.table(Table_cost, file = "Table_cost_case_1_Banking.log")

U_year = t(U)

filenamePLOT = paste("Plot_s_1_case_1_Banking",  ".pdf", sep = "")

pdf(filenamePLOT)

plot(U_year[3,], U_year[1,], type = 's', , col='gray', lwd = 1, xlab = 'Divisibility player', ylab = 'Efficiency player', ylim = c(min(U_year[1,]), max(U_year[2,])) )

polygon.step(U_year[3,], U_year[1,1:19], rep(min(U_year[1,]),19), col='gray')

lines(U_year[3,], U_year[2,], type = 's', col = "black", lty = 4, lwd = 4, xlab = 'Divisibility player', ylab = 'Efficiency player')

lines(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 2 , col = 2)

points(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 6 , col = 2)

dev.off()




###################################################################
###################################################################
#
# Compare DEA vs FDH vs NDA
#
###################################################################
###################################################################

setwd("Banking_Data")

data = read.table("Aly,Grabowski,Pasurka&Rangan-1990REStat-data.txt", header = FALSE)

Names = c('ID', 'BankID', 'X_1', 'X_2', 'X_3', 'P_1', 'P_2', 'P_3', 'Y_1', 'Y_2', 'Y_3', 'Y_4', 'Y_5', 'branch-dummy')

length(Names)

names(data) = Names

XX = data[,3:5]
YY = data[,9:13]
P_I = data[,6:8] 

m = ncol(YY)
n = ncol(XX)
k = nrow(XX)

kk = 10

k = nrow(XX)
m = ncol(YY)
n = ncol(XX)	

M_max = 20
d_lev = 2
sss = 1

aa = 0.3

DIR_M = c(rep(">=", n + m), "=", rep("<=", k))
DIR = c(rep(">=", q_mat(sss)), "<=")
DIR_M_relax = c(DIR_M, DIR_M, rep("<=", k), rep(">=", k) )
INTEGER_var = (k + q_mat(1-sss) + 1):(k + q_mat(1-sss) + k)

LL_second_block0 = cbind(diag(k), matrix(0,k,q_mat(1-sss)) )
LL_second_block1 = cbind(-diag(k), matrix(0,k,q_mat(1-sss)) )

LL_second_block = cbind(LL_second_block0, LL_second_block1)

d_val = as.matrix(c(1, seq(0.90, 1.1, length = d_lev)))
M_val = as.matrix(1:M_max)


setwd("Results/Banking")

#for (uuu in 1:k){
#	print(Efficiency(ddd = 1.1, MMM = 2, kk =uuu, ss = 1))
#}

Pi_DEA = list()
Pi_FDH = list()
Pi_NDA = list()
pb = list()

for(uuu in 1:k)
{
	#-------------------------------
	print(uuu)
	#-------------------------------

	Pi_DEA[[uuu]] = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = 1000, kk = uuu, ss = sss)
	Pi_NDA[[uuu]] = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = 2, kk = uuu, ss = sss)
	Pi_FDH[[uuu]] = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = 1, kk = uuu, ss = sss)
	pb[[uuu]] = -as.numeric(p_mat(sss,uuu)%*%b_mat(1-sss,uuu))
}

PPi_DEA =rep(0, k)
PPi_FDH =rep(0, k)
PPi_NDA =rep(0, k)

for(uuu in 1:k)
{
	PPi_DEA[uuu] = min(Pi_DEA[[uuu]]) 
	PPi_FDH[uuu] = min(Pi_FDH[[uuu]]) 
	PPi_NDA[uuu] = min(Pi_NDA[[uuu]]) 

}

summary(cbind(PPi_DEA, PPi_FDH, PPi_NDA, PPi_NDA/PPi_FDH, PPi_NDA/PPi_DEA))


A = matrix(0, 3, 3)

for(uuu in 1:k)
{

	#-------------------------------
	print(uuu)
	#-------------------------------

	if( (min(Pi_DEA[[uuu]]) < pb[[uuu]])  )
	{
		A[1,1] = A[1,1] + 1
	}

	if( (min(Pi_DEA[[uuu]]) < pb[[uuu]]) && (min(Pi_NDA[[uuu]]) < pb[[uuu]]) )
	{
		A[1,2] = A[1,2] + 1
	}

	if( (min(Pi_DEA[[uuu]]) < pb[[uuu]]) && (min(Pi_FDH[[uuu]]) < pb[[uuu]]) )
	{
		A[1,3] = A[1,3] + 1
	}


	#-------------------------------

	if( (min(Pi_NDA[[uuu]]) < pb[[uuu]]) && (min(Pi_DEA[[uuu]]) < pb[[uuu]]) )
	{
		A[2,1] = A[2,1] + 1
	}

	if( (min(Pi_NDA[[uuu]]) < pb[[uuu]])  )
	{
		A[2,2] = A[2,2] + 1
	}

	if( (min(Pi_NDA[[uuu]]) < pb[[uuu]]) && (min(Pi_FDH[[uuu]]) < pb[[uuu]]) )
	{
		A[2,3] = A[2,3] + 1
	}


	#-------------------------------


	if( (min(Pi_FDH[[uuu]]) < pb[[uuu]]) && (min(Pi_DEA[[uuu]]) < pb[[uuu]]) )
	{
		A[3,1] = A[3,1] + 1
	}



	if( (min(Pi_FDH[[uuu]]) < pb[[uuu]]) && (min(Pi_NDA[[uuu]]) < pb[[uuu]]) )
	{
		A[3,2] = A[3,2] + 1
	}


	if( (min(Pi_FDH[[uuu]]) < pb[[uuu]])  )
	{
		A[3,3] = A[3,3] + 1
	}


}

save(A, file = "Banking_Efficiency_comparaison.RData")



