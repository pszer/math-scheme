;;; real imag pair representation of complex numbers
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
