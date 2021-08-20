:- module rmath.

:- interface.
:- import_module float, int.
% 1/pi
:- func m_1_pi = float.
% 1/sqrt(2pi)
:- func m_1_sqrt_2pi = float.
% 2*pi
:- func m_2pi = float.
% 2/pi
:- func m_2_pi = float.
% 2/sqrt(pi)
:- func m_2_sqrtpi = float.
% e
:- func m_e = float.
% ln(10)
:- func m_ln10 = float.
% ln(2)
:- func m_ln2 = float.
% log(2*pi)
:- func m_ln_2pi = float.
% log(2*pi)/2
:- func m_ln_sqrt_2pi = float.
% log(pi)/2
:- func m_ln_sqrt_pi = float.
% log(pi/2)/2
:- func m_ln_sqrt_pid2 = float.
% log10(e)
:- func m_log10e = float.
% log10(2)
:- func m_log10_2 = float.
% log2(e)
:- func m_log2e = float.
% pi
:- func m_pi = float.
% pi/2
:- func m_pi_2 = float.
% pi/4
:- func m_pi_4 = float.
% 1/sqrt(2)
:- func m_sqrt1_2 = float.
% sqrt(2)
:- func m_sqrt2 = float.
% sqrt(2/pi)
:- func m_sqrt_2dpi = float.
% sqrt(3)
:- func m_sqrt_3 = float.
% sqrt(32)
:- func m_sqrt_32 = float.
% sqrt(pi)
:- func m_sqrt_pi = float.
% Exponential random numbers
:- func exp_rand = float.
% Normal random numbers
:- func norm_rand = float.
% Uniform random numbers
:- func unif_rand = float.
% R_pow function: r_pow(X,Y)
:- func r_pow(float,float) = float.
% R_pow_di function: r_pow_di(X,Y)
:- func r_pow_di(float,int) = float.
% bessel_i: bessel_i(X,Nu,Scaled)
:- func bessel_i(float,float,float) = float.
% bessel_j: bessel_j(X,Nu)
:- func bessel_j(float,float) = float.
% bessel_k: bessel_k(X,Nu,Scaled)
:- func bessel_k(float,float,float) = float.
% bessel_y: bessel_y(X,Nu)
:- func bessel_y(float,float) = float.
% beta: beta(X,Y)
:- func beta(float,float) = float.
% choose: choose(N,K)
:- func choose(float,float) = float.
% cospi: cospi(X)
:- func cospi(float) = float.
% Beta probability density function: dbeta(P,Shape1,Shape2,Log)
:- func dbeta(float,float,float,int) = float.
% Binomial probability density function: dbinom(P,Size,Prob,Log)
:- func dbinom(float,float,float,int) = float.
% Cauchy probability density function: dcauchy(P,Location,Scale,Log)
:- func dcauchy(float,float,float,int) = float.
% Chi-squared probability density function: dchisq(P,Df,Log)
:- func dchisq(float,float,int) = float.
% Exponential probability density function: dexp(P,Rate,Log)
:- func dexp(float,float,int) = float.
% F probability density function: df(P,Df1,Df2,Log)
:- func df(float,float,float,int) = float.
% Gamma probability density function: dgamma(P,Shape,Scale,Log)
:- func dgamma(float,float,float,int) = float.
% Geometric probability density function: dgeom(P,Prob,Log)
:- func dgeom(float,float,int) = float.
% Hypergeometric probability density function: dhyper(P,M,N,K,Log)
:- func dhyper(float,float,float,float,int) = float.
% digamma: digamma(X)
:- func digamma(float) = float.
% Log-normal probability density function: dlnorm(P,Meanlog,Sdlog,Log)
:- func dlnorm(float,float,float,int) = float.
% Logistic probability density function: dlogis(P,Location,Scale,Log)
:- func dlogis(float,float,float,int) = float.
% Non-central beta probability density function: dnbeta(X,Shape1,Shape2,Ncp,Log)
:- func dnbeta(float,float,float,float,int) = float.
% Negative Binomial probability density function: dnbinom(P,Size,Prob,Log)
:- func dnbinom(float,float,float,int) = float.
% Non-central chi-squared probability density function: dnchisq(P,Df,Ncp,Log)
:- func dnchisq(float,float,float,int) = float.
% Non-central F probability density function: dnf(X,Df1,Df2,Ncp,Log)
:- func dnf(float,float,float,float,int) = float.
% Normal probability density function: dnorm(P,Mean,Sd,Log)
:- func dnorm(float,float,float,int) = float.
% Non-central Student t probability density function: dnt(X,Df,Ncp,Log)
:- func dnt(float,float,float,int) = float.
% Poisson probability density function: dpois(P,Lambda,Log)
:- func dpois(float,float,int) = float.
% Wilcoxon signed rank probability density function: dsignrank(P,N,Log)
:- func dsignrank(float,float,int) = float.
% T probability density function: dt(P,Df,Log)
:- func dt(float,float,int) = float.
% Uniform probability density function: dunif(P,Min,Max,Log)
:- func dunif(float,float,float,int) = float.
% Weibull probability density function: dweibull(P,Shape,Scale,Log)
:- func dweibull(float,float,float,int) = float.
% Wilcoxon rank sum probability density function: dwilcox(P,M,N,Log)
:- func dwilcox(float,float,float,int) = float.
% fmax2: fmax2(X,Y)
:- func fmax2(float,float) = float.
% fmin2: fmin2(X,Y)
:- func fmin2(float,float) = float.
% fprec: fprec(X,Y)
:- func fprec(float,float) = float.
% fround: fround(X,Y)
:- func fround(float,float) = float.
% fsign: fsign(X,Y)
:- func fsign(float,float) = float.
% ftrunc: ftrunc(X)
:- func ftrunc(float) = float.
% gammafn: gammafn(X)
:- func gammafn(float) = float.
% get_seed: get_seed(A,B)
:- impure pred get_seed(uint32::out,uint32::out) is semidet.
% imax2: imax2(X,Y)
:- func imax2(int,int) = int.
% imin2: imin2(X,Y)
:- func imin2(int,int) = int.
% lbeta: lbeta(X,Y)
:- func lbeta(float,float) = float.
% lchoose: lchoose(N,K)
:- func lchoose(float,float) = float.
% Accurate log(gamma(x+1)) for small x (0 < x < 0.5): lgamma1p(X)
:- func lgamma1p(float) = float.
% lgammafn: lgammafn(X)
:- func lgammafn(float) = float.
% log(1 + exp(x)): log1pexp(X)
:- func log1pexp(float) = float.
% Accurate log(1+x) - x (care for small x): log1pmx(X)
:- func log1pmx(float) = float.
% log (exp (logx) + exp (logy)): logspace_add(Logx,Logy)
:- func logspace_add(float,float) = float.
% log (exp (logx) - exp (logy)): logspace_sub(Logx,Logy)
:- func logspace_sub(float,float) = float.
% Beta cumulative density function: pbeta(Q,Shape1,Shape2,Lower,Log)
:- func pbeta(float,float,float,int,int) = float.
% Binomial cumulative density function: pbinom(Q,Size,Prob,Lower,Log)
:- func pbinom(float,float,float,int,int) = float.
% Cauchy cumulative density function: pcauchy(Q,Location,Scale,Lower,Log)
:- func pcauchy(float,float,float,int,int) = float.
% Chi-squared cumulative density function: pchisq(Q,Df,Lower,Log)
:- func pchisq(float,float,int,int) = float.
% pentagamma: pentagamma(X)
:- func pentagamma(float) = float.
% Exponential cumulative density function: pexp(Q,Rate,Lower,Log)
:- func pexp(float,float,int,int) = float.
% F cumulative density function: pf(Q,Df1,Df2,Lower,Log)
:- func pf(float,float,float,int,int) = float.
% Gamma cumulative density function: pgamma(Q,Shape,Scale,Lower,Log)
:- func pgamma(float,float,float,int,int) = float.
% Geometric cumulative density function: pgeom(Q,Prob,Lower,Log)
:- func pgeom(float,float,int,int) = float.
% Hypergeometric cumulative density function: phyper(Q,M,N,K,Lower,Log)
:- func phyper(float,float,float,float,int,int) = float.
% Log-normal cumulative density function: plnorm(Q,Meanlog,Sdlog,Lower,Log)
:- func plnorm(float,float,float,int,int) = float.
% Logistic cumulative density function: plogis(Q,Location,Scale,Lower,Log)
:- func plogis(float,float,float,int,int) = float.
% Non-central beta cumulative distribution function: pnbeta(Q,Shape1,Shape2,Ncp,Lower,Log)
:- func pnbeta(float,float,float,float,int,int) = float.
% Negative Binomial cumulative density function: pnbinom(Q,Size,Prob,Lower,Log)
:- func pnbinom(float,float,float,int,int) = float.
% Non-central chi-squared cumulative density function: pnchisq(Q,Df,Ncp,Lower,Log)
:- func pnchisq(float,float,float,int,int) = float.
% Non-central F cumulative distribution function: pnf(Q,Df1,Df2,Ncp,Lower,Log)
:- func pnf(float,float,float,float,int,int) = float.
% Normal cumulative density function: pnorm(Q,Mean,Sd,Lower,Log)
:- func pnorm(float,float,float,int,int) = float.
% Non-central Student t cumulative distribution function: pnt(Q,Df,Ncp,Lower,Log)
:- func pnt(float,float,float,int,int) = float.
% Poisson cumulative density function: ppois(Q,Lambda,Lower,Log)
:- func ppois(float,float,int,int) = float.
% psigamma: psigamma(X,Y)
:- func psigamma(float,float) = float.
% Wilcoxon signed rank cumulative density function: psignrank(Q,N,Lower,Log)
:- func psignrank(float,float,int,int) = float.
% T cumulative density function: pt(Q,Df,Lower,Log)
:- func pt(float,float,int,int) = float.
% Studentised range cumulative distribution function: ptukey(Q,Nmeans,Df,Nranges,Lower,Log)
:- func ptukey(float,float,float,float,int,int) = float.
% Uniform cumulative density function: punif(Q,Min,Max,Lower,Log)
:- func punif(float,float,float,int,int) = float.
% Weibull cumulative density function: pweibull(Q,Shape,Scale,Lower,Log)
:- func pweibull(float,float,float,int,int) = float.
% Wilcoxon rank sum cumulative density function: pwilcox(Q,M,N,Lower,Log)
:- func pwilcox(float,float,float,int,int) = float.
% Beta quantile function: qbeta(P,Shape1,Shape2,Lower,Log)
:- func qbeta(float,float,float,int,int) = float.
% Binomial quantile function: qbinom(P,Size,Prob,Lower,Log)
:- func qbinom(float,float,float,int,int) = float.
% Cauchy quantile function: qcauchy(P,Location,Scale,Lower,Log)
:- func qcauchy(float,float,float,int,int) = float.
% Chi-squared quantile function: qchisq(P,Df,Lower,Log)
:- func qchisq(float,float,int,int) = float.
% Exponential quantile function: qexp(P,Rate,Lower,Log)
:- func qexp(float,float,int,int) = float.
% F quantile function: qf(P,Df1,Df2,Lower,Log)
:- func qf(float,float,float,int,int) = float.
% Gamma quantile function: qgamma(P,Shape,Scale,Lower,Log)
:- func qgamma(float,float,float,int,int) = float.
% Geometric quantile function: qgeom(P,Prob,Lower,Log)
:- func qgeom(float,float,int,int) = float.
% Hypergeometric quantile function: qhyper(P,M,N,K,Lower,Log)
:- func qhyper(float,float,float,float,int,int) = float.
% Log-normal quantile function: qlnorm(P,Meanlog,Sdlog,Lower,Log)
:- func qlnorm(float,float,float,int,int) = float.
% Logistic quantile function: qlogis(P,Location,Scale,Lower,Log)
:- func qlogis(float,float,float,int,int) = float.
% Non-central beta quantile function: qnbeta(P,Shape1,Shape2,Ncp,Lower,Log)
:- func qnbeta(float,float,float,float,int,int) = float.
% Negative Binomial quantile function: qnbinom(P,Size,Prob,Lower,Log)
:- func qnbinom(float,float,float,int,int) = float.
% Non-central chi-squared quantile function: qnchisq(P,Df,Ncp,Lower,Log)
:- func qnchisq(float,float,float,int,int) = float.
% Non-central F quantile function: qnf(P,Df1,Df2,Ncp,Lower,Log)
:- func qnf(float,float,float,float,int,int) = float.
% Normal quantile function: qnorm(P,Mean,Sd,Lower,Log)
:- func qnorm(float,float,float,int,int) = float.
% Non-central Student t quantile function: qnt(P,Df,Ncp,Lower,Log)
:- func qnt(float,float,float,int,int) = float.
% Poisson quantile function: qpois(P,Lambda,Lower,Log)
:- func qpois(float,float,int,int) = float.
% Wilcoxon signed rank quantile function: qsignrank(P,N,Lower,Log)
:- func qsignrank(float,float,int,int) = float.
% T quantile function: qt(P,Df,Lower,Log)
:- func qt(float,float,int,int) = float.
% Studentised range quantile function: qtukey(P,Nmeans,Df,Nranges,Lower,Log)
:- func qtukey(float,float,float,float,int,int) = float.
% Uniform quantile function: qunif(P,Min,Max,Lower,Log)
:- func qunif(float,float,float,int,int) = float.
% Weibull quantile function: qweibull(P,Shape,Scale,Lower,Log)
:- func qweibull(float,float,float,int,int) = float.
% Wilcoxon rank sum quantile function: qwilcox(P,M,N,Lower,Log)
:- func qwilcox(float,float,float,int,int) = float.
% Beta random numbers: rbeta(Shape1,Shape2)
:- impure func rbeta(float,float) = float.
% Binomial random numbers: rbinom(Size,Prob)
:- impure func rbinom(float,float) = float.
% Cauchy random numbers: rcauchy(Location,Scale)
:- impure func rcauchy(float,float) = float.
% Chi-squared random numbers: rchisq(Df)
:- impure func rchisq(float) = float.
% Exponential random numbers: rexp(Rate)
:- impure func rexp(float) = float.
% F random numbers: rf(Df1,Df2)
:- impure func rf(float,float) = float.
% Gamma random numbers: rgamma(Shape,Scale)
:- impure func rgamma(float,float) = float.
% Geometric random numbers: rgeom(Prob)
:- impure func rgeom(float) = float.
% Hypergeometric random numbers: rhyper(M,N,K)
:- impure func rhyper(float,float,float) = float.
% Log-normal random numbers: rlnorm(Meanlog,Sdlog)
:- impure func rlnorm(float,float) = float.
% Logistic random numbers: rlogis(Location,Scale)
:- impure func rlogis(float,float) = float.
% Negative Binomial random numbers: rnbinom(Size,Prob)
:- impure func rnbinom(float,float) = float.
% Non-central chi-squared random numbers: rnchisq(Df,Ncp)
:- impure func rnchisq(float,float) = float.
% Normal random numbers: rnorm(Mean,Sd)
:- impure func rnorm(float,float) = float.
% Poisson random numbers: rpois(Lambda)
:- impure func rpois(float) = float.
% Wilcoxon signed rank random numbers: rsignrank(N)
:- impure func rsignrank(float) = float.
% T random numbers: rt(Df)
:- impure func rt(float) = float.
% Uniform random numbers: runif(Min,Max)
:- impure func runif(float,float) = float.
% Weibull random numbers: rweibull(Shape,Scale)
:- impure func rweibull(float,float) = float.
% Wilcoxon rank sum random numbers: rwilcox(M,N)
:- impure func rwilcox(float,float) = float.
% set_seed: set_seed(A,B)
:- impure pred set_seed(uint32::in,uint32::in) is semidet.
% sign: sign(X)
:- func sign(float) = float.
% sinpi: sinpi(X)
:- func sinpi(float) = float.
% tanpi: tanpi(X)
:- func tanpi(float) = float.
% tetragamma: tetragamma(X)
:- func tetragamma(float) = float.
% trigamma: trigamma(X)
:- func trigamma(float) = float.

:- implementation.
:- pragma foreign_decl("C", "
#define MATHLIB_STANDALONE
#include \"Rmath.h\"
#define pnorm pnorm5
#define qnorm qnorm5
#define dnorm dnorm4
").
:- pragma foreign_proc("C",
    m_1_pi = (M_1_PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_1_PI_ = M_1_PI; ").
:- pragma foreign_proc("C",
    m_1_sqrt_2pi = (M_1_SQRT_2PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_1_SQRT_2PI_ = M_1_SQRT_2PI; ").
:- pragma foreign_proc("C",
    m_2pi = (M_2PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_2PI_ = M_2PI; ").
:- pragma foreign_proc("C",
    m_2_pi = (M_2_PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_2_PI_ = M_2_PI; ").
:- pragma foreign_proc("C",
    m_2_sqrtpi = (M_2_SQRTPI_::out),
    [promise_pure, will_not_call_mercury],
    "M_2_SQRTPI_ = M_2_SQRTPI; ").
:- pragma foreign_proc("C",
    m_e = (M_E_::out),
    [promise_pure, will_not_call_mercury],
    "M_E_ = M_E; ").
:- pragma foreign_proc("C",
    m_ln10 = (M_LN10_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN10_ = M_LN10; ").
:- pragma foreign_proc("C",
    m_ln2 = (M_LN2_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN2_ = M_LN2; ").
:- pragma foreign_proc("C",
    m_ln_2pi = (M_LN_2PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN_2PI_ = M_LN_2PI; ").
:- pragma foreign_proc("C",
    m_ln_sqrt_2pi = (M_LN_SQRT_2PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN_SQRT_2PI_ = M_LN_SQRT_2PI; ").
:- pragma foreign_proc("C",
    m_ln_sqrt_pi = (M_LN_SQRT_PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN_SQRT_PI_ = M_LN_SQRT_PI; ").
:- pragma foreign_proc("C",
    m_ln_sqrt_pid2 = (M_LN_SQRT_PId2_::out),
    [promise_pure, will_not_call_mercury],
    "M_LN_SQRT_PId2_ = M_LN_SQRT_PId2; ").
:- pragma foreign_proc("C",
    m_log10e = (M_LOG10E_::out),
    [promise_pure, will_not_call_mercury],
    "M_LOG10E_ = M_LOG10E; ").
:- pragma foreign_proc("C",
    m_log10_2 = (M_LOG10_2_::out),
    [promise_pure, will_not_call_mercury],
    "M_LOG10_2_ = M_LOG10_2; ").
:- pragma foreign_proc("C",
    m_log2e = (M_LOG2E_::out),
    [promise_pure, will_not_call_mercury],
    "M_LOG2E_ = M_LOG2E; ").
:- pragma foreign_proc("C",
    m_pi = (M_PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_PI_ = M_PI; ").
:- pragma foreign_proc("C",
    m_pi_2 = (M_PI_2_::out),
    [promise_pure, will_not_call_mercury],
    "M_PI_2_ = M_PI_2; ").
:- pragma foreign_proc("C",
    m_pi_4 = (M_PI_4_::out),
    [promise_pure, will_not_call_mercury],
    "M_PI_4_ = M_PI_4; ").
:- pragma foreign_proc("C",
    m_sqrt1_2 = (M_SQRT1_2_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT1_2_ = M_SQRT1_2; ").
:- pragma foreign_proc("C",
    m_sqrt2 = (M_SQRT2_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT2_ = M_SQRT2; ").
:- pragma foreign_proc("C",
    m_sqrt_2dpi = (M_SQRT_2dPI_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT_2dPI_ = M_SQRT_2dPI; ").
:- pragma foreign_proc("C",
    m_sqrt_3 = (M_SQRT_3_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT_3_ = M_SQRT_3; ").
:- pragma foreign_proc("C",
    m_sqrt_32 = (M_SQRT_32_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT_32_ = M_SQRT_32; ").
:- pragma foreign_proc("C",
    m_sqrt_pi = (M_SQRT_PI_::out),
    [promise_pure, will_not_call_mercury],
    "M_SQRT_PI_ = M_SQRT_PI; ").
:- pragma foreign_proc("C",
    exp_rand = (Exp_rand_::out),
    [promise_pure, will_not_call_mercury],
    "Exp_rand_ = exp_rand(); ").
:- pragma foreign_proc("C",
    norm_rand = (Norm_rand_::out),
    [promise_pure, will_not_call_mercury],
    "Norm_rand_ = norm_rand(); ").
:- pragma foreign_proc("C",
    unif_rand = (Unif_rand_::out),
    [promise_pure, will_not_call_mercury],
    "Unif_rand_ = unif_rand(); ").
:- pragma foreign_proc("C",
    r_pow(X::in,Y::in) = (R_pow_::out),
    [promise_pure, will_not_call_mercury],
    "R_pow_ = R_pow(X,Y);").
:- pragma foreign_proc("C",
    r_pow_di(X::in,Y::in) = (R_pow_di_::out),
    [promise_pure, will_not_call_mercury],
    "R_pow_di_ = R_pow_di(X,Y);").
:- pragma foreign_proc("C",
    bessel_i(X::in,Nu::in,Scaled::in) = (Bessel_i_::out),
    [promise_pure, will_not_call_mercury],
    "Bessel_i_ = bessel_i(X,Nu,Scaled);").
:- pragma foreign_proc("C",
    bessel_j(X::in,Nu::in) = (Bessel_j_::out),
    [promise_pure, will_not_call_mercury],
    "Bessel_j_ = bessel_j(X,Nu);").
:- pragma foreign_proc("C",
    bessel_k(X::in,Nu::in,Scaled::in) = (Bessel_k_::out),
    [promise_pure, will_not_call_mercury],
    "Bessel_k_ = bessel_k(X,Nu,Scaled);").
:- pragma foreign_proc("C",
    bessel_y(X::in,Nu::in) = (Bessel_y_::out),
    [promise_pure, will_not_call_mercury],
    "Bessel_y_ = bessel_y(X,Nu);").
:- pragma foreign_proc("C",
    beta(X::in,Y::in) = (Beta_::out),
    [promise_pure, will_not_call_mercury],
    "Beta_ = beta(X,Y);").
:- pragma foreign_proc("C",
    choose(N::in,K::in) = (Choose_::out),
    [promise_pure, will_not_call_mercury],
    "Choose_ = choose(N,K);").
:- pragma foreign_proc("C",
    cospi(X::in) = (Cospi_::out),
    [promise_pure, will_not_call_mercury],
    "Cospi_ = cospi(X);").
:- pragma foreign_proc("C",
    dbeta(P::in,Shape1::in,Shape2::in,Log::in) = (Dbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Dbeta_ = dbeta(P,Shape1,Shape2,Log);").
:- pragma foreign_proc("C",
    dbinom(P::in,Size::in,Prob::in,Log::in) = (Dbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Dbinom_ = dbinom(P,Size,Prob,Log);").
:- pragma foreign_proc("C",
    dcauchy(P::in,Location::in,Scale::in,Log::in) = (Dcauchy_::out),
    [promise_pure, will_not_call_mercury],
    "Dcauchy_ = dcauchy(P,Location,Scale,Log);").
:- pragma foreign_proc("C",
    dchisq(P::in,Df::in,Log::in) = (Dchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Dchisq_ = dchisq(P,Df,Log);").
:- pragma foreign_proc("C",
    dexp(P::in,Rate::in,Log::in) = (Dexp_::out),
    [promise_pure, will_not_call_mercury],
    "Dexp_ = dexp(P,Rate,Log);").
:- pragma foreign_proc("C",
    df(P::in,Df1::in,Df2::in,Log::in) = (Df_::out),
    [promise_pure, will_not_call_mercury],
    "Df_ = df(P,Df1,Df2,Log);").
:- pragma foreign_proc("C",
    dgamma(P::in,Shape::in,Scale::in,Log::in) = (Dgamma_::out),
    [promise_pure, will_not_call_mercury],
    "Dgamma_ = dgamma(P,Shape,Scale,Log);").
:- pragma foreign_proc("C",
    dgeom(P::in,Prob::in,Log::in) = (Dgeom_::out),
    [promise_pure, will_not_call_mercury],
    "Dgeom_ = dgeom(P,Prob,Log);").
:- pragma foreign_proc("C",
    dhyper(P::in,M::in,N::in,K::in,Log::in) = (Dhyper_::out),
    [promise_pure, will_not_call_mercury],
    "Dhyper_ = dhyper(P,M,N,K,Log);").
:- pragma foreign_proc("C",
    digamma(X::in) = (Digamma_::out),
    [promise_pure, will_not_call_mercury],
    "Digamma_ = digamma(X);").
:- pragma foreign_proc("C",
    dlnorm(P::in,Meanlog::in,Sdlog::in,Log::in) = (Dlnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Dlnorm_ = dlnorm(P,Meanlog,Sdlog,Log);").
:- pragma foreign_proc("C",
    dlogis(P::in,Location::in,Scale::in,Log::in) = (Dlogis_::out),
    [promise_pure, will_not_call_mercury],
    "Dlogis_ = dlogis(P,Location,Scale,Log);").
:- pragma foreign_proc("C",
    dnbeta(X::in,Shape1::in,Shape2::in,Ncp::in,Log::in) = (Dnbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Dnbeta_ = dnbeta(X,Shape1,Shape2,Ncp,Log);").
:- pragma foreign_proc("C",
    dnbinom(P::in,Size::in,Prob::in,Log::in) = (Dnbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Dnbinom_ = dnbinom(P,Size,Prob,Log);").
:- pragma foreign_proc("C",
    dnchisq(P::in,Df::in,Ncp::in,Log::in) = (Dnchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Dnchisq_ = dnchisq(P,Df,Ncp,Log);").
:- pragma foreign_proc("C",
    dnf(X::in,Df1::in,Df2::in,Ncp::in,Log::in) = (Dnf_::out),
    [promise_pure, will_not_call_mercury],
    "Dnf_ = dnf(X,Df1,Df2,Ncp,Log);").
:- pragma foreign_proc("C",
    dnorm(P::in,Mean::in,Sd::in,Log::in) = (Dnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Dnorm_ = dnorm(P,Mean,Sd,Log);").
:- pragma foreign_proc("C",
    dnt(X::in,Df::in,Ncp::in,Log::in) = (Dnt_::out),
    [promise_pure, will_not_call_mercury],
    "Dnt_ = dnt(X,Df,Ncp,Log);").
:- pragma foreign_proc("C",
    dpois(P::in,Lambda::in,Log::in) = (Dpois_::out),
    [promise_pure, will_not_call_mercury],
    "Dpois_ = dpois(P,Lambda,Log);").
:- pragma foreign_proc("C",
    dsignrank(P::in,N::in,Log::in) = (Dsignrank_::out),
    [promise_pure, will_not_call_mercury],
    "Dsignrank_ = dsignrank(P,N,Log);").
:- pragma foreign_proc("C",
    dt(P::in,Df::in,Log::in) = (Dt_::out),
    [promise_pure, will_not_call_mercury],
    "Dt_ = dt(P,Df,Log);").
:- pragma foreign_proc("C",
    dunif(P::in,Min::in,Max::in,Log::in) = (Dunif_::out),
    [promise_pure, will_not_call_mercury],
    "Dunif_ = dunif(P,Min,Max,Log);").
:- pragma foreign_proc("C",
    dweibull(P::in,Shape::in,Scale::in,Log::in) = (Dweibull_::out),
    [promise_pure, will_not_call_mercury],
    "Dweibull_ = dweibull(P,Shape,Scale,Log);").
:- pragma foreign_proc("C",
    dwilcox(P::in,M::in,N::in,Log::in) = (Dwilcox_::out),
    [promise_pure, will_not_call_mercury],
    "Dwilcox_ = dwilcox(P,M,N,Log);").
:- pragma foreign_proc("C",
    fmax2(X::in,Y::in) = (Fmax2_::out),
    [promise_pure, will_not_call_mercury],
    "Fmax2_ = fmax2(X,Y);").
:- pragma foreign_proc("C",
    fmin2(X::in,Y::in) = (Fmin2_::out),
    [promise_pure, will_not_call_mercury],
    "Fmin2_ = fmin2(X,Y);").
:- pragma foreign_proc("C",
    fprec(X::in,Y::in) = (Fprec_::out),
    [promise_pure, will_not_call_mercury],
    "Fprec_ = fprec(X,Y);").
:- pragma foreign_proc("C",
    fround(X::in,Y::in) = (Fround_::out),
    [promise_pure, will_not_call_mercury],
    "Fround_ = fround(X,Y);").
:- pragma foreign_proc("C",
    fsign(X::in,Y::in) = (Fsign_::out),
    [promise_pure, will_not_call_mercury],
    "Fsign_ = fsign(X,Y);").
:- pragma foreign_proc("C",
    ftrunc(X::in) = (Ftrunc_::out),
    [promise_pure, will_not_call_mercury],
    "Ftrunc_ = ftrunc(X);").
:- pragma foreign_proc("C",
    get_seed(A::out,B::out),
    [will_not_call_mercury],
    "get_seed(&A,&B);").
:- pragma foreign_proc("C",
    gammafn(X::in) = (Gammafn_::out),
    [promise_pure, will_not_call_mercury],
    "Gammafn_ = gammafn(X);").
:- pragma foreign_proc("C",
    imax2(X::in,Y::in) = (Imax2_::out),
    [promise_pure, will_not_call_mercury],
    "Imax2_ = imax2(X,Y);").
:- pragma foreign_proc("C",
    imin2(X::in,Y::in) = (Imin2_::out),
    [promise_pure, will_not_call_mercury],
    "Imin2_ = imin2(X,Y);").
:- pragma foreign_proc("C",
    lbeta(X::in,Y::in) = (Lbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Lbeta_ = lbeta(X,Y);").
:- pragma foreign_proc("C",
    lchoose(N::in,K::in) = (Lchoose_::out),
    [promise_pure, will_not_call_mercury],
    "Lchoose_ = lchoose(N,K);").
:- pragma foreign_proc("C",
    lgamma1p(X::in) = (Lgamma1p_::out),
    [promise_pure, will_not_call_mercury],
    "Lgamma1p_ = lgamma1p(X);").
:- pragma foreign_proc("C",
    lgammafn(X::in) = (Lgammafn_::out),
    [promise_pure, will_not_call_mercury],
    "Lgammafn_ = lgammafn(X);").
:- pragma foreign_proc("C",
    log1pexp(X::in) = (Log1pexp_::out),
    [promise_pure, will_not_call_mercury],
    "Log1pexp_ = log1pexp(X);").
:- pragma foreign_proc("C",
    log1pmx(X::in) = (Log1pmx_::out),
    [promise_pure, will_not_call_mercury],
    "Log1pmx_ = log1pmx(X);").
:- pragma foreign_proc("C",
    logspace_add(Logx::in,Logy::in) = (Logspace_add_::out),
    [promise_pure, will_not_call_mercury],
    "Logspace_add_ = logspace_add(Logx,Logy);").
:- pragma foreign_proc("C",
    logspace_sub(Logx::in,Logy::in) = (Logspace_sub_::out),
    [promise_pure, will_not_call_mercury],
    "Logspace_sub_ = logspace_sub(Logx,Logy);").
:- pragma foreign_proc("C",
    pbeta(Q::in,Shape1::in,Shape2::in,Lower::in,Log::in) = (Pbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Pbeta_ = pbeta(Q,Shape1,Shape2,Lower,Log);").
:- pragma foreign_proc("C",
    pbinom(Q::in,Size::in,Prob::in,Lower::in,Log::in) = (Pbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Pbinom_ = pbinom(Q,Size,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    pcauchy(Q::in,Location::in,Scale::in,Lower::in,Log::in) = (Pcauchy_::out),
    [promise_pure, will_not_call_mercury],
    "Pcauchy_ = pcauchy(Q,Location,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    pchisq(Q::in,Df::in,Lower::in,Log::in) = (Pchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Pchisq_ = pchisq(Q,Df,Lower,Log);").
:- pragma foreign_proc("C",
    pentagamma(X::in) = (Pentagamma_::out),
    [promise_pure, will_not_call_mercury],
    "Pentagamma_ = pentagamma(X);").
:- pragma foreign_proc("C",
    pexp(Q::in,Rate::in,Lower::in,Log::in) = (Pexp_::out),
    [promise_pure, will_not_call_mercury],
    "Pexp_ = pexp(Q,Rate,Lower,Log);").
:- pragma foreign_proc("C",
    pf(Q::in,Df1::in,Df2::in,Lower::in,Log::in) = (Pf_::out),
    [promise_pure, will_not_call_mercury],
    "Pf_ = pf(Q,Df1,Df2,Lower,Log);").
:- pragma foreign_proc("C",
    pgamma(Q::in,Shape::in,Scale::in,Lower::in,Log::in) = (Pgamma_::out),
    [promise_pure, will_not_call_mercury],
    "Pgamma_ = pgamma(Q,Shape,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    pgeom(Q::in,Prob::in,Lower::in,Log::in) = (Pgeom_::out),
    [promise_pure, will_not_call_mercury],
    "Pgeom_ = pgeom(Q,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    phyper(Q::in,M::in,N::in,K::in,Lower::in,Log::in) = (Phyper_::out),
    [promise_pure, will_not_call_mercury],
    "Phyper_ = phyper(Q,M,N,K,Lower,Log);").
:- pragma foreign_proc("C",
    plnorm(Q::in,Meanlog::in,Sdlog::in,Lower::in,Log::in) = (Plnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Plnorm_ = plnorm(Q,Meanlog,Sdlog,Lower,Log);").
:- pragma foreign_proc("C",
    plogis(Q::in,Location::in,Scale::in,Lower::in,Log::in) = (Plogis_::out),
    [promise_pure, will_not_call_mercury],
    "Plogis_ = plogis(Q,Location,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    pnbeta(Q::in,Shape1::in,Shape2::in,Ncp::in,Lower::in,Log::in) = (Pnbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Pnbeta_ = pnbeta(Q,Shape1,Shape2,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    pnbinom(Q::in,Size::in,Prob::in,Lower::in,Log::in) = (Pnbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Pnbinom_ = pnbinom(Q,Size,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    pnchisq(Q::in,Df::in,Ncp::in,Lower::in,Log::in) = (Pnchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Pnchisq_ = pnchisq(Q,Df,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    pnf(Q::in,Df1::in,Df2::in,Ncp::in,Lower::in,Log::in) = (Pnf_::out),
    [promise_pure, will_not_call_mercury],
    "Pnf_ = pnf(Q,Df1,Df2,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    pnorm(Q::in,Mean::in,Sd::in,Lower::in,Log::in) = (Pnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Pnorm_ = pnorm(Q,Mean,Sd,Lower,Log);").
:- pragma foreign_proc("C",
    pnt(Q::in,Df::in,Ncp::in,Lower::in,Log::in) = (Pnt_::out),
    [promise_pure, will_not_call_mercury],
    "Pnt_ = pnt(Q,Df,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    ppois(Q::in,Lambda::in,Lower::in,Log::in) = (Ppois_::out),
    [promise_pure, will_not_call_mercury],
    "Ppois_ = ppois(Q,Lambda,Lower,Log);").
:- pragma foreign_proc("C",
    psigamma(X::in,Y::in) = (Psigamma_::out),
    [promise_pure, will_not_call_mercury],
    "Psigamma_ = psigamma(X,Y);").
:- pragma foreign_proc("C",
    psignrank(Q::in,N::in,Lower::in,Log::in) = (Psignrank_::out),
    [promise_pure, will_not_call_mercury],
    "Psignrank_ = psignrank(Q,N,Lower,Log);").
:- pragma foreign_proc("C",
    pt(Q::in,Df::in,Lower::in,Log::in) = (Pt_::out),
    [promise_pure, will_not_call_mercury],
    "Pt_ = pt(Q,Df,Lower,Log);").
:- pragma foreign_proc("C",
    ptukey(Q::in,Nmeans::in,Df::in,Nranges::in,Lower::in,Log::in) = (Ptukey_::out),
    [promise_pure, will_not_call_mercury],
    "Ptukey_ = ptukey(Q,Nmeans,Df,Nranges,Lower,Log);").
:- pragma foreign_proc("C",
    punif(Q::in,Min::in,Max::in,Lower::in,Log::in) = (Punif_::out),
    [promise_pure, will_not_call_mercury],
    "Punif_ = punif(Q,Min,Max,Lower,Log);").
:- pragma foreign_proc("C",
    pweibull(Q::in,Shape::in,Scale::in,Lower::in,Log::in) = (Pweibull_::out),
    [promise_pure, will_not_call_mercury],
    "Pweibull_ = pweibull(Q,Shape,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    pwilcox(Q::in,M::in,N::in,Lower::in,Log::in) = (Pwilcox_::out),
    [promise_pure, will_not_call_mercury],
    "Pwilcox_ = pwilcox(Q,M,N,Lower,Log);").
:- pragma foreign_proc("C",
    qbeta(P::in,Shape1::in,Shape2::in,Lower::in,Log::in) = (Qbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Qbeta_ = qbeta(P,Shape1,Shape2,Lower,Log);").
:- pragma foreign_proc("C",
    qbinom(P::in,Size::in,Prob::in,Lower::in,Log::in) = (Qbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Qbinom_ = qbinom(P,Size,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    qcauchy(P::in,Location::in,Scale::in,Lower::in,Log::in) = (Qcauchy_::out),
    [promise_pure, will_not_call_mercury],
    "Qcauchy_ = qcauchy(P,Location,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    qchisq(P::in,Df::in,Lower::in,Log::in) = (Qchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Qchisq_ = qchisq(P,Df,Lower,Log);").
:- pragma foreign_proc("C",
    qexp(P::in,Rate::in,Lower::in,Log::in) = (Qexp_::out),
    [promise_pure, will_not_call_mercury],
    "Qexp_ = qexp(P,Rate,Lower,Log);").
:- pragma foreign_proc("C",
    qf(P::in,Df1::in,Df2::in,Lower::in,Log::in) = (Qf_::out),
    [promise_pure, will_not_call_mercury],
    "Qf_ = qf(P,Df1,Df2,Lower,Log);").
:- pragma foreign_proc("C",
    qgamma(P::in,Shape::in,Scale::in,Lower::in,Log::in) = (Qgamma_::out),
    [promise_pure, will_not_call_mercury],
    "Qgamma_ = qgamma(P,Shape,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    qgeom(P::in,Prob::in,Lower::in,Log::in) = (Qgeom_::out),
    [promise_pure, will_not_call_mercury],
    "Qgeom_ = qgeom(P,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    qhyper(P::in,M::in,N::in,K::in,Lower::in,Log::in) = (Qhyper_::out),
    [promise_pure, will_not_call_mercury],
    "Qhyper_ = qhyper(P,M,N,K,Lower,Log);").
:- pragma foreign_proc("C",
    qlnorm(P::in,Meanlog::in,Sdlog::in,Lower::in,Log::in) = (Qlnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Qlnorm_ = qlnorm(P,Meanlog,Sdlog,Lower,Log);").
:- pragma foreign_proc("C",
    qlogis(P::in,Location::in,Scale::in,Lower::in,Log::in) = (Qlogis_::out),
    [promise_pure, will_not_call_mercury],
    "Qlogis_ = qlogis(P,Location,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    qnbeta(P::in,Shape1::in,Shape2::in,Ncp::in,Lower::in,Log::in) = (Qnbeta_::out),
    [promise_pure, will_not_call_mercury],
    "Qnbeta_ = qnbeta(P,Shape1,Shape2,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    qnbinom(P::in,Size::in,Prob::in,Lower::in,Log::in) = (Qnbinom_::out),
    [promise_pure, will_not_call_mercury],
    "Qnbinom_ = qnbinom(P,Size,Prob,Lower,Log);").
:- pragma foreign_proc("C",
    qnchisq(P::in,Df::in,Ncp::in,Lower::in,Log::in) = (Qnchisq_::out),
    [promise_pure, will_not_call_mercury],
    "Qnchisq_ = qnchisq(P,Df,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    qnf(P::in,Df1::in,Df2::in,Ncp::in,Lower::in,Log::in) = (Qnf_::out),
    [promise_pure, will_not_call_mercury],
    "Qnf_ = qnf(P,Df1,Df2,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    qnorm(P::in,Mean::in,Sd::in,Lower::in,Log::in) = (Qnorm_::out),
    [promise_pure, will_not_call_mercury],
    "Qnorm_ = qnorm(P,Mean,Sd,Lower,Log);").
:- pragma foreign_proc("C",
    qnt(P::in,Df::in,Ncp::in,Lower::in,Log::in) = (Qnt_::out),
    [promise_pure, will_not_call_mercury],
    "Qnt_ = qnt(P,Df,Ncp,Lower,Log);").
:- pragma foreign_proc("C",
    qpois(P::in,Lambda::in,Lower::in,Log::in) = (Qpois_::out),
    [promise_pure, will_not_call_mercury],
    "Qpois_ = qpois(P,Lambda,Lower,Log);").
:- pragma foreign_proc("C",
    qsignrank(P::in,N::in,Lower::in,Log::in) = (Qsignrank_::out),
    [promise_pure, will_not_call_mercury],
    "Qsignrank_ = qsignrank(P,N,Lower,Log);").
:- pragma foreign_proc("C",
    qt(P::in,Df::in,Lower::in,Log::in) = (Qt_::out),
    [promise_pure, will_not_call_mercury],
    "Qt_ = qt(P,Df,Lower,Log);").
:- pragma foreign_proc("C",
    qtukey(P::in,Nmeans::in,Df::in,Nranges::in,Lower::in,Log::in) = (Qtukey_::out),
    [promise_pure, will_not_call_mercury],
    "Qtukey_ = qtukey(P,Nmeans,Df,Nranges,Lower,Log);").
:- pragma foreign_proc("C",
    qunif(P::in,Min::in,Max::in,Lower::in,Log::in) = (Qunif_::out),
    [promise_pure, will_not_call_mercury],
    "Qunif_ = qunif(P,Min,Max,Lower,Log);").
:- pragma foreign_proc("C",
    qweibull(P::in,Shape::in,Scale::in,Lower::in,Log::in) = (Qweibull_::out),
    [promise_pure, will_not_call_mercury],
    "Qweibull_ = qweibull(P,Shape,Scale,Lower,Log);").
:- pragma foreign_proc("C",
    qwilcox(P::in,M::in,N::in,Lower::in,Log::in) = (Qwilcox_::out),
    [promise_pure, will_not_call_mercury],
    "Qwilcox_ = qwilcox(P,M,N,Lower,Log);").
:- pragma foreign_proc("C",
    rbeta(Shape1::in,Shape2::in) = (Rbeta_::out),
    [will_not_call_mercury],
    "Rbeta_ = rbeta(Shape1,Shape2);").
:- pragma foreign_proc("C",
    rbinom(Size::in,Prob::in) = (Rbinom_::out),
    [will_not_call_mercury],
    "Rbinom_ = rbinom(Size,Prob);").
:- pragma foreign_proc("C",
    rcauchy(Location::in,Scale::in) = (Rcauchy_::out),
    [will_not_call_mercury],
    "Rcauchy_ = rcauchy(Location,Scale);").
:- pragma foreign_proc("C",
    rchisq(Df::in) = (Rchisq_::out),
    [will_not_call_mercury],
    "Rchisq_ = rchisq(Df);").
:- pragma foreign_proc("C",
    rexp(Rate::in) = (Rexp_::out),
    [will_not_call_mercury],
    "Rexp_ = rexp(Rate);").
:- pragma foreign_proc("C",
    rf(Df1::in,Df2::in) = (Rf_::out),
    [will_not_call_mercury],
    "Rf_ = rf(Df1,Df2);").
:- pragma foreign_proc("C",
    rgamma(Shape::in,Scale::in) = (Rgamma_::out),
    [will_not_call_mercury],
    "Rgamma_ = rgamma(Shape,Scale);").
:- pragma foreign_proc("C",
    rgeom(Prob::in) = (Rgeom_::out),
    [will_not_call_mercury],
    "Rgeom_ = rgeom(Prob);").
:- pragma foreign_proc("C",
    rhyper(M::in,N::in,K::in) = (Rhyper_::out),
    [will_not_call_mercury],
    "Rhyper_ = rhyper(M,N,K);").
:- pragma foreign_proc("C",
    rlnorm(Meanlog::in,Sdlog::in) = (Rlnorm_::out),
    [will_not_call_mercury],
    "Rlnorm_ = rlnorm(Meanlog,Sdlog);").
:- pragma foreign_proc("C",
    rlogis(Location::in,Scale::in) = (Rlogis_::out),
    [will_not_call_mercury],
    "Rlogis_ = rlogis(Location,Scale);").
:- pragma foreign_proc("C",
    rnbinom(Size::in,Prob::in) = (Rnbinom_::out),
    [will_not_call_mercury],
    "Rnbinom_ = rnbinom(Size,Prob);").
:- pragma foreign_proc("C",
    rnchisq(Df::in,Ncp::in) = (Rnchisq_::out),
    [will_not_call_mercury],
    "Rnchisq_ = rnchisq(Df,Ncp);").
:- pragma foreign_proc("C",
    rnorm(Mean::in,Sd::in) = (Rnorm_::out),
    [will_not_call_mercury],
    "Rnorm_ = rnorm(Mean,Sd);").
:- pragma foreign_proc("C",
    rpois(Lambda::in) = (Rpois_::out),
    [will_not_call_mercury],
    "Rpois_ = rpois(Lambda);").
:- pragma foreign_proc("C",
    rsignrank(N::in) = (Rsignrank_::out),
    [will_not_call_mercury],
    "Rsignrank_ = rsignrank(N);").
:- pragma foreign_proc("C",
    rt(Df::in) = (Rt_::out),
    [will_not_call_mercury],
    "Rt_ = rt(Df);").
:- pragma foreign_proc("C",
    runif(Min::in,Max::in) = (Runif_::out),
    [will_not_call_mercury],
    "Runif_ = runif(Min,Max);").
:- pragma foreign_proc("C",
    rweibull(Shape::in,Scale::in) = (Rweibull_::out),
    [will_not_call_mercury],
    "Rweibull_ = rweibull(Shape,Scale);").
:- pragma foreign_proc("C",
    rwilcox(M::in,N::in) = (Rwilcox_::out),
    [will_not_call_mercury],
    "Rwilcox_ = rwilcox(M,N);").
:- pragma foreign_proc("C",
    set_seed(A::in,B::in),
    [will_not_call_mercury],
    "set_seed(A,B);").
:- pragma foreign_proc("C",
    sign(X::in) = (Sign_::out),
    [promise_pure, will_not_call_mercury],
    "Sign_ = sign(X);").
:- pragma foreign_proc("C",
    sinpi(X::in) = (Sinpi_::out),
    [promise_pure, will_not_call_mercury],
    "Sinpi_ = sinpi(X);").
:- pragma foreign_proc("C",
    tanpi(X::in) = (Tanpi_::out),
    [promise_pure, will_not_call_mercury],
    "Tanpi_ = tanpi(X);").
:- pragma foreign_proc("C",
    tetragamma(X::in) = (Tetragamma_::out),
    [promise_pure, will_not_call_mercury],
    "Tetragamma_ = tetragamma(X);").
:- pragma foreign_proc("C",
    trigamma(X::in) = (Trigamma_::out),
    [promise_pure, will_not_call_mercury],
    "Trigamma_ = trigamma(X);").
