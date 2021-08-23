:- module test_rmath.

:- interface.
:- import_module io.
:- pred main(io::di, io::uo) is det.

:- implementation.
:- import_module int, float, rmath, pair, bool.

%% %% example wrapper for a random number function
%% :- impure pred wrapped_runif(float::in, float::in, float::out, io::di, io::uo).
%% wrapped_runif(Lower, Upper, U, !IO) :- U = runif(Lower,Upper).

:- type alternative ---> two_sided ; less ; greater.
:- pred poisson_ci(float::in, float::in, alternative::in, pair(float)::out).
poisson_ci(X, Conflevel, Alternative, Interval) :-
    Alpha = (1.0-Conflevel)/2.0,
    Pl = (func(Xi,Alphai) = (if Xi=0.0 then 0.0 else rmath.qgamma(Alphai,Xi, 1.0, 1, 0))),
    Pu = (func(Xi,Alphai) = rmath.qgamma(1.0-Alphai, Xi+1.0, 1.0, 1, 0)),
    Interval = (Alternative = less -> (0.0 - Pu(X, 1.0-Conflevel))
	       ;
	       Alternative = greater -> (Pl(X, 1.0-Conflevel) - 1.0)
	       ;
	       %% two_sided
	       (Pl(X,Alpha) - Pu(X, Alpha))).

:- func for_loop(func(int,int) = int, int, int, int) = int.
for_loop(Fun, I, Finish, Agg) = (if I>Finish then Agg else for_loop(Fun, I+1, Finish, Fun(I,Agg))).
:- func count(func(int) = bool, int, int) = int.
count(Predicate, Start, Finish) = Result :-
    Result = for_loop(func(I, Y) = (if Predicate(I)=yes then Y+1 else Y), Start, Finish, 0).

:- func loop1(int,float,float) = int.
loop1(Ni,M,D) = (if rmath.dpois(float(Ni),M,0)>D then loop1(Ni*2,M,D) else Ni).
:- func poisson_test(float,float,float,alternative) = float.
poisson_test(X, T, R, Alternative) = Result :-
    M = R*T,
    (Alternative = less -> Result = rmath.ppois(X,M,1,0)
    ;
    Alternative = greater -> Result = rmath.ppois(X-1.0,M,0,0)
    ;
    %% Alternative = two_sided
    (M = 0.0 -> Result = (X=0.0 -> 1.0; 0.0)
	   ;
	   (Relerr = 1.00000001,
            D = rmath.dpois(X,M,0),
	    Dstar = D * Relerr,
	    Pred = (func(I) = (if rmath.dpois(float(I),M,0) =< Dstar then yes else no)),
            (X=M -> Result = 1.0
	     ;
	     X<M ->
	     (N = loop1(ceiling_to_int(2.0*M-X),M,D),
	      Y = count(Pred, ceiling_to_int(M), N),
	      Result = rmath.ppois(X,M,1,0) + rmath.ppois(float(N)-float(Y),M,0,0))
	     ;
	     %% X>M
	     (Y = count(Pred,0,floor_to_int(M)),
	      Result = rmath.ppois(float(Y)-1.0,M,1,0) + rmath.ppois(X-1.0, M,0,0)))))).

main(!IO) :-
    runif(0.0, 1.0, U, !IO),
    runif(0.0, 1.0, U2, !IO),
    poisson_ci(5.0, 0.95, two_sided, Interval),
    P = poisson_test(5.0, 1.0, 1.0, two_sided),
    P2 = poisson_test(5.0, 10.0, 1.0, two_sided),
    io.write_line({U, U2,
		   m_e,
		   rmath.pnorm(1.96, 0.0, 1.0, 1, 0),
		   rmath.qnorm(0.975, 0.0, 1.0, 1, 0),
		   Interval,
		   P, P2
		  },
		  !IO).
