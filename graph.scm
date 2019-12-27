;;; Simple graphing tool

(import (prefix sdl2 "sdl2:"))

(define (to-int num) (inexact->exact (floor num)))

(define (graph func x y zoom w h)
	(sdl2:set-main-ready!)
	(sdl2:init! '(video))
	(define window   (sdl2:create-window! "GRAPH" 0 0 w h))
	(define renderer (sdl2:create-renderer! window))
	(define step (* 2 (/ zoom w)))

	; only look for quit event
	(sdl2:event-state-set! 'quit #t)

	(define (y-graph-to-screen gy)
		(+ (/ h 2) (/ (- y gy) step)))

	(define (draw)
		(define (graph-line screen-x graph-x last-y)
			(if (> screen-x w)
				#f
				(let ((calc-y (y-graph-to-screen (func graph-x))))
					(if (and (not (nan? calc-y)) (not (nan? last-y)))
						(sdl2:render-draw-line!
							renderer
							screen-x
							(to-int last-y)
							(+ 1 screen-x)
							(to-int calc-y)) #f)
					(graph-line (+ 1 screen-x) (+ graph-x step) calc-y))))

		(sdl2:render-draw-color-set! renderer (sdl2:make-color))
		(sdl2:render-clear! renderer)
		(sdl2:render-draw-color-set! renderer (sdl2:make-color 255 255 255))
		(graph-line 0 (+ x (- zoom) step) (y-graph-to-screen (func (- x zoom))))
		(sdl2:render-present! renderer))

	(define (main)
		(if (not (sdl2:quit-requested?))
			(main)
			#f))
	(draw)
	(sdl2:render-draw-line! renderer 0 0 600 600)
	(main)
	(sdl2:quit!))
