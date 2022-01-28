;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II
;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II
;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II  EILENBERG-MOORE-SPECTRAL-SEQUENCE-II

(IN-PACKAGE #:cat)

(PROVIDE "EILENBERG-MOORE-SPECTRAL-SEQUENCE-II")

;;We compute the tensor product of the bar complex (corresponding to WF) and the tensor product (BxF).
(DEFUN BASE-HAT-U-U (fibration
			&aux (space (sorc fibration))
			(fiber (trgt fibration))
			(bar (bar fiber))
			(spac-tnsr-fiber (tnsr-prdc space fiber)))
	(declare
		(type fibration fibration)
		(type simplicial-set space)
		(type chain-complex bar spac-tnsr-fiber)
		(type simplicial-group fiber))
	(the chain-complex
		(tnsr-prdc spac-tnsr-fiber bar)))
		
		
		
		
		
;;We compute the definition of the perturbation, corresponding to the product of the bar complex.
(DEFUN BASE-HAT-LEFT-PERTURBATION-INTR (fibration &aux (fiber (trgt fibration)))
	(declare
		(type fibration fibration)
		(type simplicial-group fiber))
	(let ((aprd (aprd fiber)))
		(flet ((rslt (degr tnpr)
			(declare
				(fixnum degr)
				(type tnpr tnpr))
			(with-tnpr (degr1 tnpr1 degr2 abar2) tnpr
				(when (zerop degr2)
					(return-from rslt (zero-cmbn (1- degr))))
				(with-tnpr (degr11 gbar11 degr12 gmsm12) tnpr1
					(let ((sign (-1-expt-n degr1))
						(brgn21 (first (abar-list abar2)))
						(bar22 (make-abar :list (rest (abar-list abar2)))))
						(declare
							(fixnum sign)
							(type brgn brgn21)
							(type abar bar22))
						(with-brgn (degr21 gmsm21) brgn21
							(let ((degr22 (- degr2 degr21)))
								(declare (fixnum degr22))
							(decf degr21)
							(let ((aprd (gnrt-? aprd (+ degr12 degr21)
										(tnpr degr12 gmsm12
										degr21 gmsm21))))
								(declare (type cmbn aprd))
								(with-cmbn (degr12+21 list) aprd
								(make-cmbn
									:degr (1- degr)
									:list 
								(mapcar
									#'(lambda (term)
									(declare (type term term))
									(with-term
										(cffc gmsm) term
										(term (* sign cffc)
										(tnpr (+ degr1 degr21)
										 (tnpr degr11 gbar11
										degr12+21 gmsm)
										degr22 bar22))))
									list)))))))))))
			(the intr-mrph #'rslt))))
			

;;Function computing the perturbation of Bar^C_*(F)(Z,Z)\otimes C_*(B)\otimes C*(F), of the coproduct.
(DEFUN BASE-HAT-LEFT-PERTURBATION (fibration)
	(declare (type fibration fibration))
	(the morphism
		(let ((base-hat-u-u (base-hat-u-u fibration)))
			(declare (type chain-complex base-hat-u-u))
			(build-mrph
				:sorc base-hat-u-u
				:trgt base-hat-u-u 
				:degr -1
				:intr (base-hat-left-perturbation-intr fibration)
				:strt :gnrt
				:orgn `(base-hat-left perturbation ,fibration)))))
				
;;Function computing the twisted tensor product: returns the reduction and the perturbation of the bcc.
;;(DEFUN BASE-SZCZARBA (fibration
;;                       &aux (space (sorc fibration))
;;                       (fiber (trgt fibration))
;;                       (twisted-crts-prdc (fibration-total fibration))
;;                       (ez (ez space fiber))
;;                       (fiber-dtau-d (fibration-dtau-d fibration)))
;;  (declare
;;   (type fibration fibration)
;;   (type simplicial-set space  twisted-crts-prdc)
;;   (type simplicial-group fiber)
;;   (type reduction ez)
;;   (type morphism fiber-dtau-d))
;;  (the (values reduction morphism) 
;;    (multiple-value-bind (szczarba bottom-perturbation)
;;        (add ez fiber-dtau-d)
;;      (with-slots (tcc f g h) szczarba
;;        (setf tcc twisted-crts-prdc
;;          (slot-value f 'sorc) twisted-crts-prdc
;;          (slot-value g 'trgt) twisted-crts-prdc
;;          (slot-value h 'sorc) twisted-crts-prdc
;;          (slot-value h 'trgt) twisted-crts-prdc))
;;      (values szczarba bottom-perturbation))))
							
							
							
							
;;Function computing the twisted tensor product of the principal fibration.
;;(DEFUN BASE-TWISTED-TNSR-PRDC (fibration &aux (szczarba (fiber-szczarba fibration)))
;;	(declare
;;		(type fibration fibration)
;;		(type reduction szczarba))
;;	(the chain-complex
;;		(bcc szczarba)))
;;Lo hago con la de eilenberg I que es lo mismo.
		
		
;;Function computing the tensor product of the bar complex and the twisted tensor product
(DEFUN BASE-HAT-U-T (fibration 
			&aux (fiber (trgt fibration))
			(bar (bar fiber))
			(twisted-tnsr-prdc (fiber-twisted-tnsr-prdc fibration))
			(base-hat-u-u (base-hat-u-u fibration)))
	(declare
		(type fibration fibration)
		(type simplicial-group fiber)
		(type chain-complex bar twisted-tnsr-prdc base-hat-u-u))
	(the chain-complex
		(let ((rslt (tnsr-prdc twisted-tnsr-prdc bar)))
			(declare (type chain-complex rslt))
			(setf (slot-value rslt 'grmd) base-hat-u-u)
			rslt)))
			
;;Function computing the twisted chain complex, with the product and the twisting operator		
(DEFUN BASE-LEFT-HMEQ-HAT (fibration 
				&aux 
				(base-hat-u-t (base-hat-u-t fibration))
				(base-hat-left-perturbation (base-hat-left-perturbation fibration)))
	(declare
		(type fibration fibration)
		(type chain-complex base-hat-u-t)
		(type morphism base-hat-left-perturbation))
	(the chain-complex
		(add base-hat-u-t base-hat-left-perturbation)))
		
		
		
;;Function computing again the twisted chain complex, with the origin.
(DEFUN EILENBERG-MOORE-BICOMPLEX-II (fibration)
	(declare (type fibration fibration))
	(let ((rslt (base-left-hmeq-hat fibration)))
		(setf (slot-value rslt 'orgn) `(Eilenberg-Moore-bicomplex-ii ,fibration))
		rslt))		
		
		
;;Method to search effective homology for the eilenberg-moore-bicomplex-ii
(DEFMETHOD SEARCH-EFHM (chcm (orgn (eql 'EILENBERG-MOORE-BICOMPLEX-II)))
	(declare (type chain-complex chcm))
	(let* ((fibration (second (orgn chcm)))
		(base (sorc fibration))
		(fiber (trgt fibration))
		(base-hat-left-perturbation (base-hat-left-perturbation fibration))
		(bar (bar fiber))
		(ttp (fiber-twisted-tnsr-prdc fibration))
		(ez (ez base fiber))
		(fiber-dtau-d (fibration-dtau-d fibration)))
		(declare
			(type fibration fibration)
			(type simplicial-set base)
			(type simplicial-group fiber)
			(type morphism base-hat-left-perturbation fiber-dtau-d)
			(type chain-complex bar ttp)
			(type reduction ez))
		(setf (slot-value bar 'efhm) (bar (efhm fiber)))
		(multiple-value-bind (rdct bottom-perturbation)
			(add ez fiber-dtau-d)
			(declare (ignore rdct)
				(type morphism bottom-perturbation))
			(setf (slot-value ttp 'efhm) (add (tnsr-prdc (efhm base) (efhm fiber)) bottom-perturbation)))
		(the homotopy-equivalence
			(add
				(tnsr-prdc (efhm ttp) (efhm bar))
				base-hat-left-perturbation))))
		

