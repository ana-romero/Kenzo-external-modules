;;;  BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM
;;;  BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM
;;;  BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM BAR-SPECTRAL-SYSTEM


(IN-PACKAGE #:cat)

(PROVIDE "bar-spectral-system")


;;This function defines the generalized filtration for the bar spectral system. It recieves a degree
;;and a generator, a tensor product of E and the bar complex. Takes the degree on the bar complex
;; degr2 (the base for the serre spectral sequence) and the number of factors in the bar (for the
;;eilenberg moore spectral sequence.
(DEFVAR TNPR-BAR-GFLIN)
(setf tnpr-bar-gflin
	#'(lambda (degr gnrt)
		(declare
			(ignore degr)
			(type tnpr gnrt))
		(with-tnpr (degr1 tnpr1 degr2 abar) gnrt
			(declare (ignore degr1))
			(let ((flin1 degr2)
			      (flin2 (with-abar (l) abar
			      		(the fixnum
			      		(length l)))))
			  (list (list (list flin1 flin2)))))))

(DEFUN CHECK-TNPR-BAR-GFLIN (fltrcm n1 &rest rest)
	(let ((gf (change-chcm-to-gflcc fltrcm (dz2) tnpr-bar-gflin '(tnpr-bar-gflin))))
		(let ((n2 (if (> (length rest) 0) (first rest) (1+ n1))))
      	(dotimes (n0 (- n2 n1))
        (let ((n (+ n0 n1)))
          (dotimes (i (length (basis gf n)))
            (if (not (gfltrd-mrph-order-p (dffr gf) n (nth i (basis gf n)) '(0 0)))
                (return-from check-tnpr-bar-gflin nil)))))
t)))


(DEFUN BAR-SPECTRAL-SYSTEM (fibration)
	(let* ((k (eilenberg-moore-bicomplex-ii fibration))
	       (ecc (rbcc (efhm k))))
		(declare
		 (type chain-complex k)
		 (type chain-complex ecc))
		(progn
			(setf k (change-chcm-to-gflcc k (dz2) tnpr-bar-gflin '(tnpr-bar-gflin)))
			(setf ecc (change-chcm-to-gflcc ecc (dz2) tnpr-bar-gflin '(tnpr-bar-gflin)))
			(the spectral-system
				(build-spectral-system ecc `(bar-spectral-system ,fibration))))))
			
