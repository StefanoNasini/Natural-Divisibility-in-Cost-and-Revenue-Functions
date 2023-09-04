

rm(list=ls(all=TRUE))
gc(reset = TRUE)

library("rPref")
library("dplyr")

library("microbenchmark")
library(Rmpfr)
library(lpSolve)
library(Matrix)



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
	#LL_second_block0 = cbind(diag(k), matrix(0,k,q_mat(1-ss)) )
	#LL_second_block1 = cbind(-diag(k), matrix(0,k,q_mat(1-ss)) )

	#LL_second_block = cbind(LL_second_block0, LL_second_block1)

	LLLP_matrix = as.matrix(rbind(LLLP_matrix, LL_second_block, LL_second_block  ))
	
	#-----------------------------------------------


	if(MMM <= 20){

		bb_k = c(b_mat(ss,kk), rep(0, q_mat(1-ss)), rep(MMM, k + 1))	
		bbb_k = c(bb_k, bb_k, aa*rep(1,k), -aa*rep(1,k) )		

		out_EDM = lp(direction = "min", objective.in = ccc_k, const.mat = LLLP_matrix, const.dir = DIR_M_relax, const.rhs = bbb_k, int.vec = INTEGER_var, timeout = 2)

		if(out_EDM$status < 2){
			#out1  = as.numeric(out_EDM$objval) 

			out1 = sum(
				out_EDM$solution[INTEGER_var]* (
					c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))
				) 
				)

		}else{
			out1  = 1.0e+20
		}

	}else{

		LLP_matrix_M = LP_mat_eq(ss, 1, ddd, slack = FALSE) 

		bb_k = c(b_mat(ss,kk), 1)

		cc_k = c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))

		out_DEA = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix_M, const.dir = DIR, const.rhs = bb_k, timeout = 2)

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

		out_EDM = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix, const.dir = DIR_M, const.rhs = bb_k, int.vec = 1:k, timeout = 2)
		
		if(out_EDM$status < 2){
			out1  = as.numeric(out_EDM$objval) 
		}else{
			out1  = 1.0e+20
		}

	}else{		

		LLP_matrix_M = LP_mat_eq(ss, 1, ddd, slack = FALSE) 

		bb_k = c(b_mat(ss,kk), 1)

		cc_k = c(-(ddd)*p_mat(ss,kk)%*%((D_mat(ss)%*%A_mat(1-ss))))

		out_DEA = lp(direction = "min", objective.in = cc_k, const.mat = LLP_matrix_M, const.dir = DIR, const.rhs = bb_k, timeout = 2)

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

# Test

#Efficiency(ddd, MMM, kk, ss)
#Relaxed_Efficiency(ddd, MMM, kk, ss)

#Efficiency(1, 5, 3, sss)
#Relaxed_Efficiency(1, 5, 3, sss)



utilities = function(MM, ss, case, PPPi, d_val){

	out1 = 0
	out1_Rel = 0

	for (kkk in 1:k){

		print(kkk)

		if (case == 0){

			if(ss == 0){

				pb = as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				tilde_a_k = abs(PPPi[[kkk]])

				Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				pb_Pi = (pb - Pi)
				pb_Pi_rel = (pb - Pi_rel)

				OUT = max(1 - pmax(0, 
						pmin( 1, pb_Pi/tilde_a_k  )
					))
				OUT_rel = max(1 - pmax(0, 
						pmin(1, pb_Pi_rel/tilde_a_k, pb_Pi/tilde_a_k )
					))
			
				#print(c(OUT, OUT_rel ))

			}else{

				tilde_a_k = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				pb = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				pb_Pi = (pb - min(Pi))
				pb_Pi_rel = (pb - min(Pi_rel))

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

				Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				pb = as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				pb_Pi = (pb - min(Pi))
				pb_Pi_rel = (pb - min(Pi_rel))

				OUT = max(1 - pmax(0, pmin(1, pb_Pi/abs(tilde_a_k))))
				OUT_rel = max(1 - pmax(0, pmin(1, pb_Pi_rel/abs(tilde_a_k))))
			
				#print(c(OUT, OUT_rel ))

			}else{

				tilde_a_k = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				Pi = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = MM, kk = kkk, ss = ss)
				Pi_rel = gmp::apply(X = d_val, MARGIN = 1, FUN = Relaxed_Efficiency, MM = MM, kk = kkk, ss = ss)

				pb = -as.numeric(p_mat(ss,kkk)%*%b_mat(1-ss,kkk))

				pb_Pi = (pb - Pi)
				pb_Pi_rel = (pb - Pi_rel)

				OUT = max(1 - pmax(0, pmin(1, pb_Pi/abs(tilde_a_k))))
				OUT_rel = max(1 - pmax(0, pmin(1, pb_Pi_rel/abs(tilde_a_k))))
			
				#print(c(OUT, OUT_rel))
			}

			out1  = out1 + OUT

			out1_Rel  = out1_Rel + OUT_rel

		}
	}

	out2 = 1 - (1/MM)
	print(MM)

	gc()

	return(c(out1/k, out1_Rel/k, out2))

}




###################################################################
###################################################################
#
# USA state-level agricultural panel data
#
###################################################################
###################################################################


#---------------------------------------------
# Load data
#---------------------------------------------


setwd("AgriculturalProduction_Data/To_import")

X1 = read.csv("Crop_output.csv", sep = ';')
X2 = read.csv("Livestock_output.csv", sep = ';')
X3 = read.csv("Other_output.csv", sep = ';')

P_I_1 = read.csv("Crop_prices.csv", sep = ';')
P_I_2 = read.csv("Livestock_prices.csv", sep = ';')
P_I_3 = read.csv("Other_output_prices.csv", sep = ';')

Y1 = read.csv("Capital_input.csv", sep = ';')
Y2 = read.csv("Land_input.csv", sep = ';')
Y3 = read.csv("Labor_input.csv", sep = ';')
Y4 = read.csv("Energy_input.csv", sep = ';')
Y5 = read.csv("Chimical_input.csv", sep = ';')

P_O_1 = read.csv("Capital_prices.csv", sep = ';')
P_O_2 = read.csv("Land_prices.csv", sep = ';')
P_O_3 = read.csv("Labor_prices.csv", sep = ';')
P_O_4 = read.csv("Energy_prices.csv", sep = ';')
P_O_5 = read.csv("Chimical_prices.csv", sep = ';')

XXX = cbind(X1[,2], X2[,2], X3[,2])
YYY = cbind(Y1[,2], Y2[2], Y3[,2], Y4[,2], Y5[,2])

k = nrow(XXX)
m = ncol(YYY)
n = ncol(XXX)

rm(XXX, YYY)
gc(reset = TRUE)

M_max = 20
d_lev = 20
sss = 0

aa = 0.3

DIR_M = c(rep(">=", n + m), "=", rep("<=", k))
DIR = c(rep(">=", q_mat(sss)), "<=")
DIR_M_relax = c(DIR_M, DIR_M, rep("<=", k), rep(">=", k) )
INTEGER_var = (k + q_mat(1-sss) + 1):(k + q_mat(1-sss) + k)

LL_second_block0 = cbind(diag(k), matrix(0,k,q_mat(1-sss)) )
LL_second_block1 = cbind(-diag(k), matrix(0,k,q_mat(1-sss)) )

LL_second_block = cbind(LL_second_block0, LL_second_block1)

setwd("Results/Plots")

Table_cost = NULL

d_val = as.matrix(c(1, seq(0.80, 1.2, length = d_lev)))
M_val = as.matrix(1:M_max)

UU = list()

PPPi = list()
PPPi_1 = list()

for(year in 2:46){

	XX = cbind(X1[,year], X2[,year], X3[,year])
	YY = cbind(Y1[,year], Y2[year], Y3[,year], Y4[,year], Y5[,year])

	P_I = cbind(P_I_1[,year], P_I_2[,year], P_I_3[,year])
	P_O = cbind(P_O_1[,year], P_O_2[year], P_O_3[,year], P_O_4[,year], P_O_5[,year])

	k = nrow(XX)
	m = ncol(YY)
	n = ncol(XX)	

	for (uuu in 1:48){
		PPPi[[uuu]] = gmp::apply(X = d_val, MARGIN = 1, FUN = Efficiency, MM = 1000, kk = uuu, ss = sss)
	}

	#utilities = function(MM, ss, case, py, PPPi, d_val)
	#utilities(3, sss, 1, PPPi, d_val)

	U = gmp::apply(X = M_val, MARGIN = 1, FUN = utilities, ss = sss, case = 0, PPPi = PPPi, d_val = d_val)

	# Compute the equilibrium point
		
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

	Table_cost = rbind(Table_cost, c(M_opt, t_opt, M_Relaxed, t_Relaxed,
							U[1,M_opt], U[2,M_Relaxed], U[3,M_opt], U[3,M_Relaxed],
							max(U[3,]), max(U[1,]),
							min(U[3,]), min(U[1,])
							))

	UU[[year]] = list("U" = U, "Table" = Table_cost)

	rm(U, tt)
	gc(reset = TRUE)

	print(Table_cost)

}

write.table(Table_cost, file = "Table_cost.log")

save(UU, file = "Full_output_cost_min.RData")


for(year in 2:46){

	U_year = UU[[year]]$U

	filenamePLOT = paste("Plot_s_0", as.character(year),  ".pdf", sep = "")

	pdf(filenamePLOT)

	plot(U_year[3,], U_year[1,], type = 's', , col='gray', lwd = 1, xlab = 'Divisibility player', ylab = 'Efficiency player', ylim = c(min(U_year[1,]), max(U_year[2,])) )

	polygon.step(U_year[3,], U_year[1,1:19], rep(min(U_year[1,]),19), col='gray')

	lines(U_year[3,], U_year[2,], type = 's', col = "black", lty = 4, lwd = 4, xlab = 'Divisibility player', ylab = 'Efficiency player')

	lines(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 2 , col = 2)

	points(c(min(U_year[3,]), max(U_year[3,])), c(min(U_year[1,]), max(U_year[2,])), lwd = 6 , col = 2)

	dev.off()

}





