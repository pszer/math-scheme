;;; zeta.scm
;;;
;;; Defines functions (zeta z) (better-zeta z)
;;; Approximations for function zeta(z) = Σ[r=0,inf] {r^-z}
;;;
;;; Definition for infinite sums for some 
;;; sequence Ai ∈ C  :  (call it partial summing)
;;; Σ[r=0,inf] Ar = lim[n->inf] Sr
;;; where Sr is the partial sum of Ar (Σ[i=0,r] Ai).
;;;
;;; zeta(z) converges completely for Real(z) > 1.0 using partial
;;; summing, diverges for any Real(z) < 1.
;;;
;;; zeta can be represented in an alternating form, which we'll call eta.
;;;
;;; Assume zeta(z) converges for z, letting us use the fact that
;;; q * (Σ[r=0,inf] {Ar}) = Σ[r=0,inf] {q*Ar}
;;; for convergent series.
;;; 
;;; eta(z) = Σ[r=0,inf] {(-1)^r * r^-z}; the alternating version of zeta
;;;
;;; P =         zeta(z) = 1^-z  +  2^-z  +  3^-z  +  4^-z  +  5^-z + ... -> inf
;;; Q = (2*2^-z)zeta(z) =        2*2^-z  +         2*4^-z  +         ... -> inf
;;; then
;;; P-Q                 = 1^-z - 2^-z + 3^-z - 4^-z + 5^-z - 6^-z ... -> inf
;;;                     = eta(z)
;;; (1-2*2^-z)zeta(z)   = eta(z)
;;;
;;; Finally, an alternative definition of zeta(z) in terms of its alternating eta
;;; form is
;;; zeta(z) = eta(z) / (1-2/2^-z)
;;;
;;; This definition is useful because the eta function converges when
;;; using 'Cesaro summation'.
;;;
;;; A Cesaro sum for some sequence Ai is
;;; lim[n->inf] Σ[r=0,n] {Sr} / n
;;; i.e. the limit of the average of the partial sums.
;;;
;;; Cesaro summation is an extension of normal partial summing
;;; i.e. if P is the set of all series evaluable by partial summing
;;;      and C is the set of all series evaluable by Cesaro summing
;;;      P is a proper subset of C.
;;;      AND if S ∈ P, then the partial sum of S = Cesaro sum of S.
;;; Cesaro summation extends the domain of the operation of
;;; infinite summing.
;;;
;;; Cesaro summation can be repeated several times by replacing Sr in the
;;; previous definition with the partial Cesaro sums Cr.
;;; From this it follows that :
;;; Let set Cn be the set of series evaluable with Cesaro sums repeated 'n' times
;;; Let set Cm be the same except repeated 'm' times, then Cn is a proper subset of Cm.
;;; i.e Repeating Cesaro summation is more powerful, [Note: limit of C[inf] != universal set] 
;;;
;;; Certain divergent series can be evaluated using Cesaro summation
;;; For example, the sequence Ai = 1 : 0 ? (even? i)
;;; A = 1,   0,   1,   0,   1,   0,   1,   0,   1,   0 (...)
;;; S = 1,   1,   2,   2,   3,   3,   4,   4,   5,   5 (...)
;;; Obviously, S doesn't converge so we can't using partial summing to evaluate
;;; the series A. But using Cesaro summation (where Cn is the average of the first
;;; n'th terms).
;;; C = 1, 1/2, 2/3, 1/2, 3/5, 1/2, 4/7, 1/2, 5/9, 1/2 (...)
;;;
;;; Ci = 1/2               : (odd?  i)
;;;    = (i+2) / (2(i+1))  : (even? i)
;;; Obviously limit of Ci -> inf is 1/2
;;; So the Cesaro sum of 1+0+1+0+1+0+1+... C= 1/2
;;; [Noted as C= instead of = to signify it's not a normal sum]
;;;
;;; eta(z) converges for Real(z) > 1.0
;;; Using Cesaro summation repeated 'n' times, the domain of
;;; eta(z) changes to Real(z) > 1.0 - n
;;; Using Cesaro summation we can evaluate eta(z) at any point
;;; z ∈ C, at least theoretically. In computing it directly
;;; we are heavily limited by the accuracy of our floating point
;;; representation.
;;; 

;;; Summing functions.
(define (sum start end term)
	(define (iter r total)
		(if (> r end)
		    total
		    (iter (+ r 1) (+ total (term r)))))
	(iter start 0))
(define (sum-complex start end term)
	(define (iter r total)
		(if (> r end)
		    total
		    (iter (+ r 1) (add-complex total (term r)))))
	(iter start (complex 0 0)))
;;; 'power' is how many times the Cesaro 'transform' is repeated (>0 int)
(define (sum-cesaro start end term power)
	(define (last seq)
		(if (null? (cdr seq))
		    (car seq)
		    (last (cdr seq))))
	(define (iter i curr-seq)
		(if (= i power)
		    (last curr-seq)
		    (iter (+ i 1) (enumerate-cesaro-sum curr-seq))))
	(iter 0 (enumerate-partial-sum start end term)))

;;; Enumerate a complex sequence as a list
(define (enumerate-sequence start end term)
  	(define (iter r result)
		(if (> r end)
		    result
		    (iter (+ r 1) (cons (term r) result))))
	(reverse (iter start '())))

;;; Returns a list of partial sums
;;; i.e (enumerate-partial-sums 0 6 identity)
;;;     '(0 1 3 6 10 15 21)
(define (enumerate-partial-sum start end term)
  	(define (iter total r result)
		(if (> r end)
		    result
		    (let ((new-total (add-complex total (term r))))
		          (iter new-total (+ r 1) (cons new-total result)))))
	(reverse (iter (complex 0 0) start '())))

(define (partial-sum sequence)
	(define (iter seq total result)
		(if (null? seq)
		    result
		    (iter (cdr seq)
			  (add-complex total (car seq))
			  (cons (add-complex total (car seq)) result))))
	(reverse (iter sequence (complex 0 0) '())))

;;; Returns a list of the partial Cesaro sums
;;; i.e (enumerate-cesaro-sum '(1 0 1 0 1 0 1 0))
;;;     '(1 1/2 2/3 1/2 3/5 1/2 4/7 1/2)
;;; NOTE: takes in a sequence of complex numbers
(define (enumerate-cesaro-sum sequence)
  	(define (iter r total seq result)
		(if (null? seq)
		    result
		    (let ((new-total (add-complex total (car seq))))
		          (iter (+ r 1)
				new-total
				(cdr seq)
				(cons (mul-complex new-total (complex (/ 1.0 r) 0))
				      result)))))
	(reverse (iter 1 (complex 0 0) sequence '())))

;;; real imag tuple representation of complex
;;; numbers
(define (complex real imag)
	(cons real imag))
(define (real c) (car c))
(define (imag c) (cdr c))
(define (neg-complex c)
	(complex (- (real c)) (- (imag c))))
(define (conj-complex c) (complex (real c) (- (imag c))))
(define (display-complex c)
	(display (real c)) (display "+")
	(display (imag c)) (display "i")
	(newline))

;;; Arithmetic operators
;;; z = a + bi
;;; z,w ∈ C
;;; w = c + di
;;; z+w = (a+c)   + (b+d)i   ∈ C
;;; z-w = (a-c)   + (b-d)i   ∈ C
;;; z*w = (ac-bd) + (ad+bc)i ∈ C
(define (add-complex a b) (complex
                               (+ (real a) (real b))
                               (+ (imag a) (imag b))))
(define (sub-complex a b) (complex
                               (- (real a) (real b))
                               (- (imag a) (imag b))))
(define (mul-complex a b . rest)
	(let ((result (complex (- (* (real a) (real b)) (* (imag a) (imag b)))
		               (+ (* (real a) (imag b)) (* (imag a) (real b))))))
		(if (null? rest)
			result
			(apply mul-complex (cons result (cons (car rest) (cdr rest)))))))

;;; z/w = (a+bi)/(c+di)
;;;     = (a+bi)(c-di)/(c+di)(c-di)
;;;     = (a+bi)(c-di)/(c^2+d^2)
;;; Denominator becomes a real number
(define (div-complex n d)
	(define (conj-mul c) ; c multiplied by its conjugate
		(define (square n) (* n n))
		(+ (square (real c)) (square (imag c))))
	(mul-complex (mul-complex n (conj-complex d))
		     (complex (/ 1 (conj-mul d)) 0)))

;;; exp-complex, raises a real number to a complex power
;;; b^z = e^(ln[b]*z) = e^(ln[b]*a)*e^(ln[b]*bi)
;;; e^(ln[b]*a)  is simply a real number
;;; e^(ln[b]*bi) = cos(ln[b]*b) + sin(ln[b]*b)i
;;; The problem of real^imag now becomes (real^real)*imag.
(define (exp-complex base exponent)
	(define const-e (exp 1))
	(let ((log-base (log base)))
		(mul-complex 
			(complex (expt const-e (* log-base (real exponent)))
			          0.0)
			(complex (cos (* log-base (imag exponent)))
			         (sin (* log-base (imag exponent)))))))

(define (zeta-term z)
	(lambda (n) (exp-complex n (neg-complex z))))
(define (eta-term z)
	(lambda (n) (if (odd? n)
			(exp-complex n (neg-complex z))
			(neg-complex (exp-complex n (neg-complex z))))))


(define (eta-coeff z) ; 1-2/2^z
	(sub-complex (complex 1 0)
		     (div-complex (complex 2 0)
				  (exp-complex 2 z))))
(define zeta-terms 8000)
(define (zeta z)
	(cond
		((number? z) (zeta (complex z 0)))
		((> (real z) 1) (sum-complex 1 zeta-terms (zeta-term z)))
		(else (div-complex (sum-cesaro 1 zeta-terms (eta-term z)
		                               (+ 1 (floor (abs (real z)))))
		                   (eta-coeff z)))
	))

;;; zeta(2) as an approximation of pi.
;;;
;;; 1. sin(x)/x  = sum [r=0,inf] (x^2r)(-1)^r / (2r+1)!
;;; 2. sin(x)/x  = prod[r=1,inf] (1-x/(r*pi))(1+x/(r*pi))
;;; 3.           = prod[r=1,inf] (1-x^2/(r*pi)^2)         ; difference of squares
;;; Expanding the infinite product brackets (.3), we get that the x^2
;;; co-efficient is -(sum[r=0,inf] {1/(r*pi)^2})
;;; The x^2 co-efficient in the series (1.) is -1/6
;;; Equating co-efficients and multiplying by -pi^2, we get
;;; sum[r=0,inf] {1/r^2} = pi^2 / 6
;;;              zeta(2) = pi^2 / 6
;;; so we can use sqrt(6*zeta(2)) as an approximation of pi.
;;; zeta(2) converges relatively quickly because 1/x^2 becomes
;;; small very quickly, but doesn't converge too quickly as an
;;; approximation of pi, because in calculating pi we use a square
;;; root on zeta(2) which effectively negates the quickly converging
;;; nature of 1/x^2 .
;(define pi-approx 
;	(sqrt (* 6 (real (sum-complex 1 2500 (zeta-term (complex 2.0 0.0)))))))

;;; A better way of calculating zeta for σ < 1 is with its
;;; 'functional form'
;;; zeta(z) = 2^z pi^(z-1) sin[(s*pi)/2] Γ(1-s) zeta(1-s) 
;;; where Γ is the gamma function.

;;; Γ(z) = int[0,inf] {x^(z-1)*e^-x dx}
(define (integral f a b dx)
	(define (add-dx x) (+ x dx))
	(define (sum term a next b)
		(accumulate add-complex (complex 0 0) term a next b))
	(define (accumulate combiner null-value term a next b)
		(define (iter total r)
			(if (> r b)
			     total
			     (iter (combiner total (term r)) (next r))))
		(iter null-value a))
	(mul-complex (sum f (+ a (/ dx 2)) add-dx b) (complex dx 0)))
(define (gamma z)
	(define m-one (sub-complex z (complex 1 0)))
	(define (term x) (mul-complex (exp-complex x m-one) (complex (exp (- x)) 0)))
	(define (square x) (* x x))
	(integral term 0 (max 100 (square (real z))) 0.05))

;;; exponential form of sin
;;; sin(z) = (e^zi - e^-zi)/2i
(define (sin-complex z)
	(define e (exp 1))
	(define zi (mul-complex z (complex 0 1)))
	(div-complex (sub-complex (exp-complex e zi)
	                          (exp-complex e (neg-complex zi)))
	             (complex 0 2)))

(define (better-zeta z)
	(define pi 3.14159265359)
	(cond
	    ((number? z) (zeta (complex z 0)))
	    ((> (real z) 0) (zeta z))
	    (else
		(mul-complex
            		(exp-complex 2 z)
	    		(exp-complex pi (sub-complex z (complex 1 0)))
	    		(sin-complex (div-complex (mul-complex z (complex pi 0)) (complex 2 0)))
	    		(gamma (sub-complex (complex 1 0) z))
	    		(zeta  (sub-complex (complex 1 0) z))))))

 
