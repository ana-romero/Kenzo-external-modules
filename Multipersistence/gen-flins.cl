(IN-PACKAGE #:cat)

(PROVIDE "gen-flins")

(DEFVAR CRPR-FLIN)
(SETF CRPR-FLIN 
  #'(lambda (degr crpr)
      (declare
       (type fixnum degr)
       (type crpr crpr))
      (the fixnum
        (- degr (length (dgop-int-ext (dgop1 crpr)))))))

(DEFVAR TNPR-FLIN)
(SETF TNPR-FLIN 
  #'(lambda (degr tnpr)
      (declare 
       (ignore degr)
       (type fixnum )
       (type tnpr tnpr))
      (the fixnum
        (degr1 tnpr))))

(DEFVAR CRPR2-GFLIN)
(SETF CRPR2-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((dgop1 (dgop1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let* ((dgop11 (dgop*dgop dgop1 (dgop1 gnrt1)))
               (dgop12 (dgop*dgop dgop1 (dgop2 gnrt1)))
               (flin1 (- degr (length (dgop-int-ext dgop11))))
               (flin2 (- degr (length (dgop-int-ext dgop12)))))
          (list (list (list flin1 flin2)))))))

(DEFVAR TNPR2-GFLIN)
(SETF TNPR2-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((flin1 (degr1 (gnrt1 gnrt)))
            (flin2 (degr2 (gnrt1 gnrt))))
        (list (list (list flin1 flin2))))))

(DEFUN 2-POINTS-ADD (p1 p2)
  (declare (type list p1 p2))
  (the list
    (list 
     (+ (first p1) (first p2))
     (+ (second p1) (second p2)))))

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
            ;;(print degr-1)
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

(DEFVAR TNPR2-GFLIN3)
(SETF TNPR2-GFLIN3 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (with-tnpr (degr1 gnrt1 degr2 gnrt2) gnrt
        (declare (ignore degr1 degr2 gnrt2))
        (let ((p (degr1 gnrt1))
              (q (degr2 gnrt1))
              (rslt nil))
          (dotimes (i (1- p))
            (push (list i (- (+ p q) (+ 1 i))) rslt))
          (push (list p q) rslt)
          (list rslt)))))

(DEFVAR TNPR2-GFLIN4)
(SETF TNPR2-GFLIN4 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (with-tnpr (degr1 gnrt1 degr2 gnrt2) gnrt
        (declare (ignore degr1 degr2 gnrt2))
        (let ((p (degr1 gnrt1))
              (q (degr2 gnrt1))
              (rslt nil))
          (dotimes (i p)
            (push (list i (- (+ p q) i)) rslt))
          (push (list p q) rslt)
          (list rslt)))))

(DEFVAR CRPR2-GFLIN2)
(SETF CRPR2-GFLIN2 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((dgop1 (dgop1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let* ((dgop11 (dgop*dgop dgop1 (dgop1 gnrt1)))
               ;;(dgop12 (dgop*dgop dgop1 (dgop2 gnrt1)))
               (flin1 (- degr (length (dgop-int-ext dgop11))))
               (flin2 (- degr (length (dgop-int-ext dgop1)))))
          (list (t-downset-list (list (- flin2 flin1) flin1)))))))

#|
(DEFVAR TNPR-CRPR-GFLIN)
(SETF TNPR-CRPR-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((degr1 (degr1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let* ((flin1 (- degr1 (length (dgop-int-ext (dgop1 gnrt1)))))
              (flin2 (- degr1 flin1)))
          (list (t-downset-list (list flin2 flin1)))))))
|#

(DEFVAR TNPR-CRPR-GFLIN)
(SETF TNPR-CRPR-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((degr1 (degr1 gnrt))
            (gnrt1 (gnrt1 gnrt)))
        (let ((flin1 (- degr1 (length (dgop-int-ext (dgop1 gnrt1)))))
              (flin2 (- degr1 (length (dgop-int-ext (dgop2 gnrt1))))))
          (list (list (list flin1 flin2)))))))



(DEFVAR TNPR2-GFLIN2)
(SETF TNPR2-GFLIN2 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((flin1 (degr1 (gnrt1 gnrt)))
            (flin2 (degr2 (gnrt1 gnrt))))
        (list (t-downset-list (list flin2 flin1))))))

(DEFVAR 3TNPR-GFLINA)
(SETF 3TNPR-GFLINA #'(lambda (degr gnrt)
                      (declare 
                       (ignore degr)
                       (type tnpr gnrt))
                      (let* ((flin1 (degr2 (gnrt1 gnrt)))
                             (flin2 (first (first (funcall tnpr2-gflin2 flin1 (gnrt1 gnrt)))))
                             )
                        (list (t3-downset-list (nconc (list flin1) flin2))))))

(DEFVAR 3TNPR-GFLIN)
(SETF 3TNPR-GFLIN #'(lambda (degr gnrt)
                      (declare 
                       (ignore degr)
                       (type tnpr gnrt))
                      (let* ((flin1 (degr1 (gnrt1 (gnrt1 gnrt))))
                             (flin2 (degr2 (gnrt1 (gnrt1 gnrt))))
                             (flin3 (degr2 (gnrt1 gnrt))))
                        (list (t3-downset-list (list flin3 flin2 flin1))))))



