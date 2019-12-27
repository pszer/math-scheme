;;; derivative.scm
;;; (deriv expression var)            - Performs symbolic derivation on an expression
;;; (calc-deriv expression var point) - Calculates a deriative at a point 

(define (calc-deriv exp var point)
	(define e 2.718281828459045)
	(define (substitute l)
		(cond
			((eq? var l) point)
			((pair? l) (cons (substitute (car l)) (substitute (cdr l))))
			(else l)))
	(let ((d (substitute (deriv exp var))))
		(display d)
		(eval d user-initial-environment)))

(define (deriv exp var)
	(define (variable? x) (symbol? x))
	(define (same-variable? v1 v2)
		(and (variable? v1) (variable? v2) (eq? v1 v2)))

	(define (=number? exp num) (and (number? exp) (= exp num)))
	(define (is-quotient? exp) (and (pair? exp) (eq? (car exp) '/)))
	(define (make-sum     a1 a2)
		(cond
			((=number? a1 0) a2)
			((=number? a2 0) a1)
			((and (number? a1) (number? a2))
				(+ a1 a2))
			(else (list '+ a1 a2))))
	(define (make-product m1 m2)
		(cond
			((or (=number? m1 0) (=number? m2 0)) 0)
			((=number? m1 1) m2)
			((=number? m2 1) m1)
			((and (is-quotient? m1) (equal? m2 (caddr m1))) 1)
			((and (is-quotient? m2) (equal? m1 (caddr m2))) 1)
			((and (number? m1) (number? m2))
				(* m1 m2))
			(else (list '* m1 m2))))
	(define (make-quotient n d)
		(cond
			((=number? n 0) 0)
			((=number? d 1) n)
			((equal? n d)   1)
			((and (number? n) (number? d))
				(/ n d))
			(else (list '/ n d))))
	(define (make-exponent b e)
		(cond
			((=number? e 1) b)
			((=number? e 0) 1)
			((=number? b 0) 0)
			((=number? b 1) 1)
			(else (list 'expt b e))))
	(define (make-log arg)
		(list 'log arg))

	;;; list of trig derivative pairs. 
	(define trig-derivs
		'((sin (cos x))                  ; sin x -> cos x
		  (cos (- (sin x)))              ; cos x -> -sin x
		  (tan (expt (/ 1 (cos x)) 2))   ; tan x -> sec^2 x
		 ))
	(define (get-trig-deriv func argument)
		(define (search seq)
			(cond
			    ((null? seq) '())
			    ((eq? func (caar seq)) (cadar seq))
			    (else (search (cdr seq)))))
		(define (tree-map proc tree)
			(map (lambda (sub-tree)
			        (if (pair? sub-tree)
			            (tree-map proc sub-tree)
	        		    (proc sub-tree)))
					tree))
		(tree-map
			(lambda (p) (if (eq? p 'x) argument p))
			(search trig-derivs)))

	(define (get-exponent-deriv base power var)
		(cond ((and (eq? base var) (number? power))  ; x^R
			(make-product power
			              (make-exponent base (- power 1))))
		      ((eq? base 'e)                         ; e^...
			(make-product (deriv power var)
			              (make-exponent 'e power))
			)
		      (else (make-product (deriv (make-product (make-log base) power) var)
		                          (make-exponent base power)))))

	(define (get-log-deriv argument var)
		(make-quotient (deriv argument var) argument))

	(define (sum?      x) (and (pair? x) (eq? (car x) '+)))
	(define (product?  x) (and (pair? x) (eq? (car x) '*)))
	(define (exponent? x) (and (pair? x) (eq? (car x) 'expt)))
	(define (log?      x) (and (pair? x) (eq? (car x) 'log)))
	(define (trig? x)
		(define (test seq)
			(cond ((null? seq) #f)
			      ((eq? (car x) (caar seq)) #t)
			      (else (test (cdr seq)))))
		(if (pair? x)
			(test trig-derivs)
			#f))

	(define (addend s) (cadr  s))
	(define (augend s)
		(let ((rest (cddr s)))
			(if (null? (cdr rest))
				(car rest)
				(cons '+ rest))))

	(define (multiplier   p) (cadr  p))
	(define (multiplicand p)
		(let ((rest (cddr p)))
			(if (null? (cdr rest))
				(car rest)
				(cons '* rest))))

	(define (exp-base  e) (cadr  e))
	(define (exp-power e) (caddr e))
	(cond
		((number? exp) 0)
		((variable? exp) (if (same-variable? exp var) 1 0))
		((sum? exp) (make-sum (deriv (addend exp) var)
		                      (deriv (augend exp) var)))
		((product? exp)
		 (make-sum
			(make-product (multiplier exp)
			              (deriv (multiplicand exp) var))
			(make-product (deriv (multiplier exp) var)
			              (multiplicand exp))))
		((exponent? exp)
			;(make-product (exp-power exp)
			;              (make-exponent (exp-base exp) (- (exp-power exp) 1))))
			(get-exponent-deriv (exp-base exp) (exp-power exp) var))
		((trig? exp)
			(make-product (deriv (cadr exp) var)
			              (get-trig-deriv (car exp) (cadr exp))))
		((log? exp)
			(get-log-deriv (cadr exp) var))
		(else
			(error "unknown expression type: DERIV" exp))))
