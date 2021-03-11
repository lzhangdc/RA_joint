varDelta <- function(p2, p){
	(2*p^2*(p2-p^2+p*(1-p)) - 4*p*(1-p)*p2 + (1-p2)*p2)
}

delta_clt <- function(g, delta0){
	n <- length(g)
	p2_hat <- mean(g==2); p_hat <- mean(g)/2
	delta_hat <- p2_hat - p_hat^2
	v_hat <- varDelta(p2_hat, p_hat)/n
	zStat <- (delta_hat - delta0)/sqrt(v_hat)
	return(c(delta_hat, zStat))
}

Tbeta <- function(gr, gs){
	r <- length(gr); s <- length(gs)
	pr <- mean(gr)/2; ps <- mean(gs)/2
	p2r <- mean(gr==2); p2s <- mean(gs==2)
	delta_r <- p2r-pr^2
	delta_s <- p2s-ps^2
	T1 <- pr-ps
	T2 <- delta_r - delta_s
	V1 <- (delta_r+pr*(1-pr))/2/r + (delta_s+ps*(1-ps))/2/s
	V2 <- varDelta(p2r, pr)/r + varDelta(p2s, ps)/s
	VC <- (p2r-pr^2)*(1-2*pr)/r + (p2s-ps^2)*(1-2*ps)/s
	s_fn <- matrix(c(T1, T2), ncol=1)
	sigma <- matrix(c(V1, VC, VC, V2), ncol=2)
	#print(sigma)
	t_stat <- t(s_fn)%*%solve(sigma)%*%s_fn
	return(t_stat)
}


geno_fn <- function(gr, gs){
	yr <- rep(1, length(gr))
	ys <- rep(0, length(gs))
	g <- factor(as.character(c(gr, gs)))
	y <- c(yr, ys)
	fit1 <- glm(y~g, family='binomial')
	fit0 <- glm(y~1, family='binomial')
	anova(fit1, fit0, test='Chisq')$'Pr(>Chi)'[2]
}


HWE_test <- function(geno_counts){
        n2 <- geno_counts[3]; n1 <- geno_counts[2]
        n <- sum(geno_counts)
        p <- (n2+n1/2)/n
        p2 <- n2/n
        delta <- p2-p^2
        rho <- delta/p/(1-p)
        p_HWE <- pchisq(n*rho^2, df=1, lower.tail=F)
        return(c(p, delta, p_HWE))
}

info_mat <- function(alpha, beta, r0, r1, r2, s0, s1, s2){
	n0 <- r0+s0; n1 <- r1+s1; n2 <- r2+s2
	r <- r0+r1+r2
	f0 <- exp(alpha)/(1+exp(alpha))
	f1 <- exp(alpha+beta)/(1+exp(alpha+beta))
	f2 <- exp(alpha+2*beta)/(1+exp(alpha+2*beta))
	T1 <- -n0*f0*(1-f0)-n1*f1*(1-f1)-n2*f2*(1-f2)
	T2 <- -n1*f1*(1-f1)-2*n2*f2*(1-f2)
	T3 <- -n1*f1*(1-f1)-4*n2*f2*(1-f2)
	L1 <- r-n0*f0-n1*f1-n2*f2
	L2 <- r1 + 2*r2 -n1*f1-2*n2*f2
	# L_vec <- matrix(c(L1, L2), ncol=1)
	# j_mat <- solve(matrix(c(T1, T2, T2, T3), ncol=2))
	matrix(c(T3*L1-T2*L2, T1*L2-T2*L1)/(T1*T3-T2^2), ncol=1)
}


nw_find <- function(r0, r1, r2, s0, s1, s2, nIter=10, crit = 10^(-10)){
	para_vec <- matrix(ncol=nIter+1, nrow=2)
	para_vec[,1] <- matrix(c(0,0), ncol=1)
	for(i in 1:nIter){
		temp_crit <- info_mat(alpha=para_vec[1,i], beta=para_vec[2,i], r0, r1, r2, s0, s1, s2)
		para_vec[,i+1] <- para_vec[,i] - temp_crit
		if(abs(temp_crit[1,1]) >= crit & abs(temp_crit[1,1]) >= crit){
			para_vec[,i+1] <- para_vec[,i] - temp_crit
		}else{
			return(para_vec[,i])
		}
	}
	return(para_vec[,(nIter+1)])
}


fast_ra_joint <- function(control_g, case_g){
        s1 <- control_g[2]; s2 <- control_g[3]; r1 <- case_g[2]; r2 <- case_g[3]
        s <- sum(control_g); r <- sum(case_g)
        p2r <- r2/r; p2s <- s2/s
        pr <- (2*r2+r1)/r/2; ps <- (2*s2+s1)/s/2
        delta_r <- p2r-pr^2
        delta_s <- p2s-ps^2
        beta <- pr-ps
        delta <- delta_r - delta_s
        V1 <- (delta_r+pr*(1-pr))/2/r + (delta_s+ps*(1-ps))/2/s
        V2 <- varDelta(p2r, pr)/r + varDelta(p2s, ps)/s
        VC <- delta_r*(1-2*pr)/r + delta_s*(1-2*ps)/s
        (beta^2*V2-2*beta*delta*VC+delta^2*V1)/(V1*V2-VC^2)
}

geno_wald_fast <- function(control_g, case_g){
        s0 <- control_g[1]; s1 <- control_g[2]; s2 <- control_g[3]
        r0 <- case_g[1]; r1 <- case_g[2]; r2 <- case_g[3]
        s <- sum(control_g); r <- sum(case_g)
        p1r <- r1/r; p2r <- r2/r; p0r <- 1-p1r-p2r
        p1s <- s1/s; p2s <- s2/s; p0s <- 1-p1s-p2s
        V0 <- 1/r0 + 1/s0; V1 <- 1/r1 + 1/s1; V2 <- 1/r2 + 1/s2
        T1 <- (log(p1r/p1s)-log(p2r/p2s))^2/(V1+V2+V1*V2/V0)
        T2 <- (log(p1r/p1s)-log(p0r/p0s))^2/(V0+V1+V0*V1/V2)
        T3 <- (log(p0r/p0s)-log(p2r/p2s))^2/(V0+V2+V0*V2/V1)
        stat <- T1 + T2 + T3
        return(stat)
        }


fast_add_wald <- function(control_g, case_g){
        r0 <- case_g[1]; r1 <- case_g[2]; r2 <- case_g[3]
        s0 <- control_g[1]; s1 <- control_g[2]; s2 <- control_g[3]
        est_vec <- nw_find(r0, r1, r2, s0, s1, s2)
        alpha <- est_vec[1]; beta <- est_vec[2]
        n0 <- r0+s0; n1 <- r1+s1; n2 <- r2+s2
        f0 <- exp(alpha)/(1+exp(alpha))
        f1 <- exp(alpha+beta)/(1+exp(alpha+beta))
        f2 <- exp(alpha+2*beta)/(1+exp(alpha+2*beta))
        T1 <- n0*f0*(1-f0)+n1*f1*(1-f1)+n2*f2*(1-f2)
        T2 <- n1*f1*(1-f1)+2*n2*f2*(1-f2)
        T3 <- n1*f1*(1-f1)+4*n2*f2*(1-f2)
        var_beta <- T3-T2^2/T1
        r <- sum(r0, r1, r2); s <- sum(s0, s1, s2)
        y <- rep(c(0, 1), c(s, r))
        g <- c(rep(0:2, control_g), rep(0:2, case_g))
        beta^2*var_beta 
}


fast_add_score <- function(case_g, control_g){
	r <- sum(case_g); s <- sum(control_g); n <- r+s
	pr <- (case_g[2]+2*case_g[3])/2/r
	ps <- (control_g[2]+2*control_g[3])/2/s
	pool_g <- case_g + control_g
	p <- (pool_g[2] + pool_g[3]*2)/2/n
	delta <- pool_g[3]/n-p^2
	denom <- (p*(1-p)+delta)*(1/r+1/s)/2
	(pr-ps)^2/denom
}

geno_score_fast <- function(case_g, control_g){
	r <- sum(case_g); s <- sum(control_g); n <- r+s
	case_freq <- case_g/r
	control_freq <- control_g/s
	pool_g <- case_g + control_g
	pool_freq <- pool_g/n
	sum(sapply(1:3, function(x) (case_freq[x]-control_freq[x])^2/pool_freq[x]))*r*s/n
}


fast_MI_joint_HWE <- function(gSNP, pheno){
        rs_table <- table(gSNP, pheno)
        if(nrow(rs_table) < 3){
                #print(ncol(rs_table))
                return(rep(NA, 14))
        }
        if(any(rs_table < 5)){
                #print(rs_table)
                return(rep(NA, 14))
        }
        if(nrow(rs_table) == 3){
                controlG <- rs_table[,1]
                caseG <- rs_table[,2]
        }
        RA_joint <- fast_ra_joint(controlG, caseG); RA_joint_pval <- pchisq(RA_joint, df=2, lower.tail = F)
        geno_stat <- geno_wald_fast(controlG, caseG); geno_pval <- pchisq(geno_stat, df=2, lower.tail=F)
        geno_score_stat <- geno_score_fast(controlG, caseG); geno_score_pval <- pchisq(geno_score_stat, df=2, lower.tail=F)
	add_wald_stat <- fast_add_wald(control_g=controlG, case_g=caseG); add_pval <- pchisq(add_wald_stat, df=1, lower.tail=F)
        add_score_stat <- fast_add_score(control_g=controlG, case_g=caseG); add_score_pval <- pchisq(add_score_stat, df=1, lower.tail=F)
	control_HWE <- HWE_test(rs_table[,1]) 
        case_HWE <- HWE_test(rs_table[,2])
        joint_HWE <- HWE_test(rs_table[,1]+rs_table[,2])
        return(c(RA_joint_pval, geno_pval, geno_score_pval, add_pval, add_score_pval, control_HWE, case_HWE, joint_HWE))
}
