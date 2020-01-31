; solves the SAT problem
; given a boolean expression of ANDed clauses
;   c0 ^ c1 ^ c2 ^ ... ^ ck
; where ci = (v0 v v1 v ... v vk) vi=0,1
; is there an assignment to the variables that makes
; the expression evaluate to 1.

; variables = list of symbols
; clauses = list of lists of symbols
; ((x1 x2 x3) (¬x1 ¬x2 x3) (¬x3 ¬x2))
; = (x1 v x2 v x3) ^ (¬x1 v ¬x2 v x3) ^ (¬x3 v ¬x2)
(define (sat variables clauses)
  (let ((env (make-top-level-environment)))
	(define (define-variables! variables)
	  (unless (null? variables)
		(environment-define env (car variables) #f)
		(define-variables! (cdr variables))))
	(define (create-expression clauses)
	  (define (iter clauses result)
		(define (format-clause clause)
		  (display clause) (newline)
		  (if (null? clause)
			'()
			(let ((symbol-str (symbol->string (car clause))))
			  (if (char=? #\¬ (string-ref symbol-str 0))
				(cons (list 'not (string->symbol (string-tail symbol-str 1)))
					  (format-clause (cdr clause)))
				(cons (car clause) (format-clause (cdr clause)))))))
		(if (null? clauses)
		  result
		  (iter (cdr clauses)
				(cons (cons 'or (format-clause (car clauses))) result))))
	  (cons 'and (iter clauses '())))
	(define-variables! variables)
	(let ((expression (create-expression clauses))
		  (result '()))
	  (display expression) (newline)
	  (define (add-assignment)
		(set! result
		  (cons (map (lambda (sym) (cons sym (environment-lookup env sym)))
					 variables)
				result)))
	  (define (explore variables)
		(if (null? variables)
		  (when (eval expression env)
			(add-assignment))
		  (begin (environment-assign! env (car variables) #f)
				 (explore (cdr variables))
				 (environment-assign! env (car variables) #t)
				 (explore (cdr variables)))))
	  (explore variables)
	  result)))
