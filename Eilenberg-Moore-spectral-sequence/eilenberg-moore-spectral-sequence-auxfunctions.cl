;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS
;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS
;;;  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS  EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS

(IN-PACKAGE #:cat)

(PROVIDE "EILENBERG-MOORE-SPECTRAL-SEQUENCE-AUXFUNCTIONS")


(DEFUN ZERO-SMMR (smst1 &optional (smst2 smst1) (degr 0))
  (declare (type simplicial-set smst1 smst2)
           (type fixnum degr))
  (the simplicial-mrph
    (build-smmr
     :sorc smst1 :trgt smst2 :degr degr
     :sintr #'(lambda (dmns gmsm)
                (declare (ignore gmsm)
                         (type fixnum dmns))
                (absm (mask (- dmns degr)) (bspn smst2)))
     :orgn `(zero-smmr ,smst1 ,smst2 ,degr))))


(DEFMETHOD CMPS ((mrph1 simplicial-mrph) (mrph2 simplicial-mrph) &optional strt )
  (declare (ignore strt))
  (the simplicial-mrph
    (with-slots ((sorc1 sorc) (trgt1 trgt) (degr1 degr)) mrph1
      (declare
       (type simplicial-set sorc1 trgt1)
       (fixnum degr1))
      (with-slots ((sorc2 sorc) (trgt2 trgt) (degr2 degr)) mrph2
        (declare
         (type simplicial-set sorc2 trgt2)
         (fixnum degr2))
        (unless (eq sorc1 trgt2)
          (error "In 2MRPH-CMPS, the morphisms ~A and ~A may not be composed (cf source and target)."
            mrph1 mrph2))        
        (when (or (eq (first (orgn mrph1)) 'zero-smmr)
                  (eq (first (orgn mrph2)) 'zero-smmr))
          (return-from cmps
            (zero-smmr sorc2 trgt1 (+ degr1 degr2))))
        (when (eq (first (orgn mrph1)) 'idnt-smmr)
          (return-from cmps mrph2))
        (when (eq (first (orgn mrph2)) 'idnt-smmr)
          (return-from cmps mrph1))
        (build-smmr
         :sorc sorc2 :trgt trgt1 :degr (+ degr1 degr2)
         :sintr #'(lambda (dmns gmsm)
                    (declare
                     (type fixnum dmns)
                     (type gmsm gmsm))
                    (the absm
                      (? mrph1 (+ dmns degr2) (? mrph2 dmns gmsm))))
         :orgn `(2mrph-cmps ,mrph1 ,mrph2 ))))))


(DEFUN TWOP-INCL-INTR (twop)
  (declare (type fibration twop))
  (flet ((rslt (dmns gmsm)
               (absm 0 (crpr  (1- (expt 2 dmns)) (bspn (sorc twop)) 0 gmsm))
               ))
    (the intr-mrph #'rslt)))


(DEFUN TWOP-INCL (twop)
  (let* ((fiber (trgt twop))
         (total (fibration-total twop)))
    (build-smmr :sorc fiber :trgt total :degr 0 :sintr (twop-incl-intr twop)
                :orgn `(fibration-inclusion ,twop))))