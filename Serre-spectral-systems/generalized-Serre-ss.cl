;; GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS
;; GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS
;; GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS   GENERALIZED-SERRE-SS


;;; FILTRATIONS FOR THE DIFFERENT SPACES

(in-package :user)


(setf CRPR-FLIN 
  #'(lambda (degr crpr)
      (declare
       (type fixnum degr)
       (type crpr crpr))
      (the fixnum
        (- degr (length (dgop-int-ext (dgop1 crpr)))))))


(setf TNPR-FLIN 
  #'(lambda (degr tnpr)
      (declare 
       (ignore degr)
       (type fixnum )
       (type tnpr tnpr))
      (the fixnum
        (degr1 tnpr))))


(DEFUN T-DOWNSET-LIST (p)
  (declare (type list p))
  (let ((p1 (first p))
        (p2 (second p))
        (rslt nil))
    (progn
      (mapcar #'(lambda (i)
                  (push (2-points-add p (list i (- i))) rslt))
        (nreverse (<a-b> 1 p2)))
      (push p rslt)
      (mapcar #'(lambda (j)
                  (push (2-points-add p (list (1- (- j)) j)) rslt))
        (<a-b> 1 (1- p1)))
      rslt
      )))


(DEFUN TWOP-PROJ-1-INTR (twop)
  (declare (type fibration twop))
  (flet ((rslt (dmns gmsm)
               (absm 0 (crpr   0 gmsm (1- (expt 2 dmns)) (bspn (trgt twop))))
               ))
    (the intr-mrph #'rslt)))


(DEFUN CHECK-TWOP (twop degr)
  (let* ((sorc (sorc twop))
         (trgt (trgt twop))
         (degr-1 (1- degr)))
    (dotimes (n (length (basis sorc degr)))
      (let ((gn (nth n (basis sorc degr))))
        (progn
          ;;(print gn)
          (dotimes (j (1- degr))
            ;;(print j)
            (let ((dtg (face trgt j (1- degr) (? twop  degr gn )))
                  (tdg (tw-a-sintr3 (sintr twop) degr-1 (face sorc j degr gn) (bspn trgt))))
              (if (and (= (dgop dtg) (dgop tdg))
                       (eq :equal (cmpr trgt (gmsm dtg)  (gmsm tdg))))
                  () 
                (progn
                  (print gn)
                  (print j)
                  (print 'nok)))))
          (let* ((dd-1tg (face trgt degr-1 degr-1 (? twop  degr gn )))
                 (tddg (tw-a-sintr3 (sintr twop) degr-1 (face sorc degr degr gn) (bspn trgt)))
                 (tddg-1 (grin trgt (1- degr-1) tddg))
                 (tdd-1g (tw-a-sintr3 (sintr twop) degr-1 (face sorc degr-1 degr gn) (bspn trgt)))
                 (tddg-1+tdd-1g (grml trgt (1- degr-1) (crpr (dgop tddg-1) (gmsm tddg-1) (dgop tdd-1g) (gmsm tdd-1g)))))
            (if (and (= (dgop dd-1tg) (dgop tddg-1+tdd-1g))
                     (eq :equal (cmpr trgt (gmsm dd-1tg)  (gmsm tddg-1+tdd-1g))))
                ()
              (progn
                (print gn)
                (print n)
                (print degr-1)
                (print 'nok)))))
        ))))


(DEFUN DZN-TRANSLATE (n ds p)
  (mapcar #'(lambda (q)
              (mapcar #'(lambda (i)
                          (+ (nth i q) (nth i p)))
                (<a-b< 0 n)))
    ds))


(DEFUN GEN-FLTRD-INDEX-p (gfltrcm degr gnrt gen-fltr-index)
  (declare (type GENERALIZED-FILTERED-CHAIN-COMPLEX gfltrcm)
           (type fixnum degr ))
  (let* ((pos (pos gfltrcm))
         (rslt nil)
         (gflin (gen-flin gfltrcm degr gnrt)))
    (declare (type list gflin))
    (mapcar #'(lambda (el)
                (if (and (or (eq :equal (pocmpr pos el gen-fltr-index))
                             (eq :less (pocmpr pos el gen-fltr-index)))
                         (not (eq (first rslt) gnrt)))
                    (setf rslt 't)))
      gflin)
    rslt))


(DEFUN GFLTRD-MRPH-ORDER-P (mrph degr gnrt1 order)
  (let* ((degrm (degr mrph))
         (cmbn (? mrph degr gnrt1))
         (cmbnl (cmbn-list cmbn))
         (rslt 't)
         (sorc (sorc mrph))
         (trgt (trgt mrph))
         (gflin (gen-flin sorc degr gnrt1)))
    (mapcar #'(lambda (flin)
                (let ((flin+order (dzn-translate (length order) flin order))
                      (b 't))
                  (mapcar #'(lambda (term)
                              (if (not (GEN-FLTRD-INDEX-p trgt (+ degr degrm) (gnrt term) flin+order))
                                  (setf b nil)))
                    cmbnl)
                  (if (not b) (setf rslt nil))
                  ))
      gflin)
    rslt))


(setf TNPR-CRPR-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((degr1 (degr1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let* ((flin1 (- degr1 (length (dgop-int-ext (dgop1 gnrt1)))))
               (flin2 (- degr1 flin1)))
          (list (t-downset-list (list flin2 flin1)))))))


(setf CRPR2-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((dgop1 (dgop1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let* ((dgop11 (dgop*dgop dgop1 (dgop1 gnrt1)))
               ;;(dgop12 (dgop*dgop dgop1 (dgop2 gnrt1)))
               (flin1 (- degr (length (dgop-int-ext dgop11))))
               (flin2 (- degr (length (dgop-int-ext dgop1))))
               )
          (list (t-downset-list (list (- flin2 flin1) flin1)))))))


(setf TNPR2-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((flin1 (degr1 (gnrt1 gnrt)))
            (flin2 (degr2 (gnrt1 gnrt))))
        (list (t-downset-list (list flin2 flin1))))))


(setf TNPR3-GFLIN #'(lambda (degr gnrt)
                      (declare 
                       (ignore degr)
                       (type tnpr gnrt))
                      (let* ((flin1 (degr1 (gnrt1 (gnrt1 gnrt))))
                             (flin2 (degr2 (gnrt1 (gnrt1 gnrt))))
                             (flin3 (degr2 (gnrt1 gnrt)))
                             )
                        (list (t3-downset-list (list flin3 flin2 flin1))))))


(setf CRPR3-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let* ((dgop1 (dgop1 gnrt))
             (gnrt1 (gnrt1 gnrt))
             (dgop11 (dgop*dgop dgop1 (dgop1 gnrt1)))
             (gnrt11 (gnrt1 gnrt1))
             (dgop111 (dgop*dgop dgop11 (dgop1 gnrt11)))
             )
        (let* (
               (flin1 (- degr (length (dgop-int-ext dgop111))))
               (flin2 (- degr (length (dgop-int-ext dgop11))))
               (flin3 (- degr (length (dgop-int-ext dgop1))))
               )
          (list (t3-downset-list (list (- flin3 flin2) (- flin2 flin1) flin1)))))))



;;; FUNCTIONS FOR THE GENERALIZED SERRE SPECTRAL SEQUENCE


(DEFUN 2-POINTS-ADD (p1 p2)
  (declare (type list p1 p2))
  (the list
    (list 
     (+ (first p1) (first p2))
     (+ (second p1) (second p2)))))


(DEFUN E2-2GSPSQ-BASIS-DVS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 1 -2))))
         (b (t-downset-list (2-points-add point (list 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-basis-dvs gfltrcm z q p b degr)))


(DEFUN E2-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 1 -2))))
         (b (t-downset-list (2-points-add point (list 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN E2-2GSPSQ-GNRTS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 1 -2))))
         (b (t-downset-list (2-points-add point (list 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN E2-EFF-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 1 -2))))
         (b (t-downset-list (2-points-add point (list 0 1)))))
    (declare (type list point p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN Ee1-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 0 -1))))
         (b (t-downset-list (2-points-add point (list 1 0)))))
    (declare (type list point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN Ee1-EFF-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 0 -1))))
         (b (t-downset-list (2-points-add point (list 1 0)))))
    (declare (type list point p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN Ee1-2GSPSQ-GNRTS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 0 -1))))
         (b (t-downset-list (2-points-add point (list 1 0)))))
    (declare (type list point p q z b))
    (gen-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN Ee1-EFF-2GSPSQ-GNRTS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z (t-downset-list (2-points-add point (list 0 -1))))
         (b (t-downset-list (2-points-add point (list 1 0)))))
    (declare (type list point p q z b))
    (gen-eff-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN E1-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN E1-EFF-2GSPSQ-GROUP (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN E1-2GSPSQ-GNRTS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN E1-EFF-2GSPSQ-GNRTS (gfltrcm p1 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr))
  (let* ((point (list p1 p2))
         (p (t-downset-list point))
         (q (t-downset-list (2-points-add point (list 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-eff-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN E1-EFF-3GSPSQ-GROUP (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN E1-3GSPSQ-GROUP (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z q)
         (b p))
    (declare (type list point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN A-DOWNSET-LIST-3 (p1 p2 p3 degr)
  (let ((rslt nil))
    (progn
      (push (list p1 p2 p3) rslt)
      (if (> p2 0) (push (list p1 (1- p2) (* 2 degr)) rslt))
      (if (> p1 0) (push (list (1- p1) (* 2 degr) (* 2 degr)) rslt))
      rslt)))


(DEFUN A-DOWNSET-LIST-2 (p1 p2 degr)
  (let ((rslt nil))
    (progn
      (push (list p1 p2 (* 2 degr)) rslt)
      (if (> p1 0) (push (list (1- p1) (* 2 degr) (* 2 degr)) rslt))
      rslt)))


(DEFUN A-DOWNSET-LIST-1 (p1 degr)
  (list  (list p1 (* 2 degr) (* 2 degr))))


(DEFUN LEXCON-EFF-3GSPSQ-GROUP-3 (gfltrcm p1 p2 p3 r degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 r degr))
  (let* ((p (a-downset-list-3 p1 p2 p3 degr))
         (q (a-downset-list-3 p1 p2 (1- p3) degr))
         (z (a-downset-list-3 p1 p2 (- p3 r) degr))
         (b (a-downset-list-3 p1 p2 (1- (+ p3 r)) degr)))
    (declare (type list p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))

(DEFUN LEXCON-EFF-3GSPSQ-GROUP-2 (gfltrcm p1 p2 r degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 r degr))
  (let* ((p (a-downset-list-2 p1 p2 degr))
         (q (a-downset-list-2 p1 (1- p2) degr))
         (z (a-downset-list-2 p1  (- p2 r) degr))
         (b (a-downset-list-2 p1  (1- (+ p2 r)) degr)))
    (declare (type list p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN LEXCON-EFF-3GSPSQ-GROUP-1 (gfltrcm p1 r degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 r degr))
  (let* ((p (a-downset-list-1 p1  degr))
         (q (a-downset-list-1  (1- p1) degr))
         (z (a-downset-list-1   (- p1 r) degr))
         (b (a-downset-list-1   (1- (+ p1 r)) degr)))
    (declare (type list p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN 2-nPOINTS-ADD (p1 p2)
  (declare (type list p1 p2))
  (the list
    (mapcar #'(lambda (i)
                (+ (nth i p1) (nth i p2)))
      (<a-b< 0 (length p1)))))


(DEFUN T3-DOWNSET-LIST (p)
  (declare (type list p))
  (let* ((p1 (first p))
         (p2 (second p))
         (p3 (third p))
         (N (+ (+ p1 p2) p3))
         (rslt nil))
    (declare
     (type fixnum p1 p2 p3)
     (type list rslt))
    (the list
      (progn
        (mapcar #'(lambda (i)
                    (mapcar #'(lambda (j)
                                (push (list i j (- (- N i) j)) rslt))
                      (nreverse (<a-b> 0 (- N i)))))
          (nreverse (<a-b> (+ 1 p1) N)))
        (mapcar #'(lambda (j)
                    (push (list p1 j (- (- N p1) j)) rslt))
          (nreverse (<a-b> p2 (+ p2 p3))))
        (if (> p1 0)
            (mapcar #'(lambda (j)
                        (push (list (1- p1) j (- (+ p2 p3) j)) rslt))
              (nreverse (<a-b> 0 (1- p2)))))
        (mapcar #'(lambda (i)
                    (mapcar #'(lambda (j)
                                (push (list i j (- (- (1- N) i) j)) rslt))
                      (nreverse (<a-b> 0 (- (1- N) i)))))
          (nreverse (<a-b> 0 (- p1 2))))
        
        rslt))))


(DEFUN E2-3GSPSQ-BASIS-DVS (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add point (list 0 1 -2))))
         (b (t3-downset-list (2-npoints-add point (list 0 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-basis-dvs gfltrcm z q p b degr)))


(DEFUN E2-3GSPSQ-GROUP (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add point (list 0 1 -2))))
         (b (t3-downset-list (2-npoints-add point (list 0 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN E2-3GSPSQ-GNRTS (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add point (list 0 1 -2))))
         (b (t3-downset-list (2-npoints-add point (list 0 0 1)))))
    (declare (type list point p q z b))
    (gen-spsq-gnrts gfltrcm z q p b degr)))


(DEFUN E2-EFF-3GSPSQ-GROUP (gfltrcm p1 p2 p3 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr))
  (let* ((point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add point (list 0 1 -2))))
         (b (t3-downset-list (2-npoints-add point (list 0 0 1)))))
    (declare (type list point p q z b))
    (gen-eff-spsq-group gfltrcm z q p b degr)))


(DEFUN NEXT-E2-3GSPSQ-GROUP-TRANS (gfltrcm p1 p2 p3 point2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr)
   (type list point2))
  (let* ((mpoint2 (mapcar #'(lambda (i)
                              (- i))
                    point2))
         (point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (dzn-translate 3 (t3-downset-list (2-npoints-add point (list 0 1 -2))) mpoint2))
         (b (dzn-translate 3 (t3-downset-list (2-npoints-add point (list 0 0 1))) point2)))
    (declare (type list mpoint2 point p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN NEXT-E2-3GSPSQ-GROUP-SUM (gfltrcm p1 p2 p3 point2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr)
   (type list point2))
  (let* ((mpoint2 (mapcar #'(lambda (i)
                              (- i))
                    point2))
         (point (list p1 p2 p3))
         (p (t3-downset-list point))
         (q (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add mpoint2 (2-npoints-add point (list 0 1 -2)))))
         (b (t3-downset-list (2-npoints-add point2 (2-npoints-add point (list 0 0 1))))))
    (declare (type list point mpoint2 p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN FINAL-2GSPSQ-GROUP (gfltrcm degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   )
  (let* ((point (list degr degr))
         (p (t-downset-list point))
         (q (t-downset-list (list -1 -1)))
         (z q)
         (b p))
    (declare (type list p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN NEXT-E2-2GSPSQ-GROUP-TRANS (gfltrcm p1 p2 point2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 degr)
   (type list point2))
  (let* ((mpoint2 (mapcar #'(lambda (i)
                              (- i))
                    point2))
         (point (list p1 p2 ))
         (p (t-downset-list point))
         (s (t-downset-list (2-points-add point (list 1 -1))))
         (z (dzn-translate 2 (t-downset-list (2-points-add point (list 0 -2))) mpoint2))
         (b (dzn-translate 2 (t-downset-list (2-points-add point (list 0 1))) point2)))
    (declare (type list point mpoint2 p s z b))
    (gen-spsq-group gfltrcm z s p b degr)))


(DEFUN LEXCON-2GSPSQ-GROUP (gfltrcm p1 p2 point2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2  degr)
   (type list point2))
  (let* ((mpoint2 (mapcar #'(lambda (i)
                              (- i))
                    point2))
         (point (list p1 p2 ))
         (p (t-downset-list point))
         (s (t-downset-list (2-npoints-add point (list 1 -1))))
         (z (t-downset-list (2-npoints-add mpoint2 (2-npoints-add point (list 1 -2)))))
         (b (t-downset-list (2-npoints-add point2 (2-npoints-add point (list 0 1))))))
    (declare (type list point mpoint2 p s z b))
    (gen-spsq-group gfltrcm z s p b degr)))


(DEFUN LEXCON-3GSPSQ-GROUP (gfltrcm p1 p2 p3 point2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum p1 p2 p3 degr)
   (type list point2))
  (let* ((mpoint2 (mapcar #'(lambda (i)
                              (- i))
                    point2))
         (point (list p1 p2 p3))
         (p (t3-downset-list point))
         (s (t3-downset-list (2-npoints-add point (list 0 1 -1))))
         (z (t3-downset-list (2-npoints-add mpoint2 (2-npoints-add point (list 0 1 -2)))))
         (b (t3-downset-list (2-npoints-add point2 (2-npoints-add point (list 0 0 1))))))
    (declare (type list point mpoint2 p s z b))
    (gen-spsq-group gfltrcm z s p b degr)))


(DEFUN FINAL-3GSPSQ-GROUP (gfltrcm degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr))
  (let* ((point (list degr degr degr))
         (p (t3-downset-list point))
         (q (t3-downset-list (list -1 -1 -1)))
         (z q)
         (b p))
    (declare (type list p q z b))
    (gen-spsq-group gfltrcm z q p b degr)))


(DEFUN E2-GSPSQ-GROUP (gfltrcm l degr)
  (if (= (length l) 2) (e2-2gspsq-group gfltrcm (first l) (second l) degr)
    (e2-3gspsq-group gfltrcm (first l) (second l) (third l) degr)))


(DEFUN LEXCON-GSPSQ-GROUP (gfltrcm l l2 degr)
  (if (= (length l) 2) (lexcon-2gspsq-group gfltrcm (first l) (second l) l2 degr)
    (lexcon-3gspsq-group gfltrcm (first l) (second l) (third l) l2 degr)))

(DEFUN FINAL-GSPSQ-GROUP (gfltrcm degr)
  (let ((orgn (orgn (second (orgn (pos gfltrcm))))))
    (if (eq 'z2 (first orgn)) (final-2gspsq-group gfltrcm degr)
      (let ((n (second orgn)))
        (if (= 2 n) (final-2gspsq-group gfltrcm degr)
          (final-3gspsq-group gfltrcm degr))))))





(DEFUN E1-EFF-GSPSQ-GROUP (gfltrcm l degr)
  (if (= (length l) 2) (E1-EFF-2GSPSQ-GROUP gfltrcm (first l) (second l) degr)
    (E1-EFF-3GSPSQ-GROUP gfltrcm (first l) (second l) (third l) degr)))

