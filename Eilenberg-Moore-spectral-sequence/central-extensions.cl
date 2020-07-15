;; CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS
;; CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS
;; CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS   CENTRAL-EXTENSIONS

(IN-PACKAGE #:cat)

(PROVIDE "central-extensions")


;;;; CARTESIAN PRODUCT OF GROUPS


(DEFUN GRCRPR-P (object)
  (declare (type any object))
  (the boolean
    (and (consp object)
         (eql :grcrpr (car object))
         (consp (cdr object))
         (typep (cadr object) 'gnrt)
         (typep (cddr object) 'gnrt)
         )))


(DEFTYPE GRCRPR () '(satisfies grcrpr-p))


#|
(DEFMETHOD PRINT-KEYCONS ((car (eql :grcrpr)) cdr stream)
  (declare
   (cons cdr)
   (stream stream))
  (the (eql t)
    (progn
      (setf cdr (cons car cdr))
      (with-grcrpr (grcrpr-e1 grcrpr-e2) cdr
        (format stream
            "<Group-CrPr ~A ~A>"
          grcrpr-e1 grcrpr-e2)
        cdr)
      t)))
|#

(DEFMACRO grcrpr (g1 g2)
  `(cons :grcrpr (cons ,g1 ,g2)))

(DEFMACRO grcrpr-e1 (grcrpr)
  `(cadr ,grcrpr))

(DEFMACRO grcrpr-e2 (grcrpr)
  `(cddr ,grcrpr))

(DEFMACRO WITH-GRCRPR ((grcrpr-e1 grcrpr-e2) grcrpr . body)
  `(let (,@(if grcrpr-e1 `((,grcrpr-e1 (grcrpr-e1 ,grcrpr))) nil)
            ,@(if grcrpr-e2 `((,grcrpr-e2 (grcrpr-e2 ,grcrpr))) nil))
     (declare
      (type gnrt ,@(if grcrpr-e1 `(,grcrpr-e1) nil) ,@(if grcrpr-e2 `(,grcrpr-e2) nil)))
     ,@body))

#|
 
  (macroexpand '(grcrpr a b))
 
  (setf grcrpr (grcrpr 'a (grcrpr 'b 'c))))
  (macroexpand-1 '(with-grcrpr (gnrt1 gnrt2) grcrpr
                    (progn gnrt1 gnrt2)))
  (macroexpand-1 '(with-grcrpr (gnrt1 gnrt2) grcrpr
                    (statement-1)
                    (statement-2)))
|#

(DEFUN GR-CRTS-PRDC-ELEMENTS (basis1 basis2)
  (declare
   (type basis basis1 basis2))
  (when (or (eq basis1 :locally-effective)
            (eq basis2 :locally-effective))
    (return-from gr-crts-prdc-elements :locally-effective))
  (the list
    (mapcan #'(lambda (e1)
                (mapcar #'(lambda (e2)
                            (grcrpr e1 e2))
                  basis2))
      basis1)))


(DEFUN GR-CRTS-PRDC-CMPR (cmpr1 cmpr2)
  (declare (type cmprf cmpr1 cmpr2))
  (flet ((rslt (grcrpr1 grcrpr2)
               (declare (type grcrpr grcrpr1 grcrpr2))
               (the cmpr
                 (let ((e11 (cadr grcrpr1))
                       (e21 (cadr grcrpr2))
                       (e12 (cddr grcrpr1))
                       (e22 (cddr grcrpr2)))
                   
                   (lexico
                    (funcall cmpr1 e11 e21)
                    (funcall cmpr2 e12 e22))))))
    (the cmprf #'rslt)))


(DEFUN GR-CRTS-PRDC-MULT (mult1 mult2)
  (declare 
   (type function mult1 mult2))
  (flet ((rslt (grcrpr1 grcrpr2)
               (with-grcrpr (e11 e12) grcrpr1
                 (with-grcrpr (e21 e22) grcrpr2
                   (grcrpr (funcall mult1 e11 e21) (funcall mult2 e12 e22))))))
    #'rslt))


(DEFUN GR-CRTS-PRDC-INV (inv1 inv2)
  (declare 
   (type function inv1 inv2))
  (flet ((rslt (grcrpr)
               (with-grcrpr (e1 e2) grcrpr
                 (grcrpr (funcall inv1 e1) (funcall inv2 e2 )))))
    #'rslt))


(DEFUN GR-CRTS-PRDC (group1 group2)
  (declare
   (type group group1 group2))
  (with-slots ((elements1 elements) (cmpr1 cmpr) (mult1 mult) (inv1 inv) (nullel1 nullel) (orgn1 orgn)) group1
    (with-slots ((elements2 elements) (cmpr2 cmpr) (mult2 mult) (inv2 inv) (nullel2 nullel) (orgn2 orgn)) group2
      (let ((elements (GR-CRTS-PRDC-ELEMENTS elements1 elements2))
            (cmpr (gr-crts-prdc-cmpr cmpr1 cmpr2))
            (mult (gr-crts-prdc-mult mult1 mult2))
            (inv (gr-crts-prdc-inv inv1 inv2))
            (nullel (grcrpr nullel1 nullel2))
            (orgn `(gr-crts-prdc ,group1 ,group2)))
        (let ((rslt (build-group :elements elements :cmpr cmpr :mult mult :inv inv :nullel nullel :orgn orgn)))
          (if (and (typep group1 'ab-group) (typep group2 'ab-group))
              (change-class rslt 'ab-group)))))))


#|
 (setf gr (cyclicgroup 5))
 (setf crpr (gr-crts-prdc gr gr))
 (setf grcrpr1 (grcrpr 3 4) grcrpr2 (grcrpr 1 2))
 (? crpr grcrpr1 grcrpr2)
 (funcall (inv1 crpr) grcrpr1)
|#

;;; EFFECTIVE HOMOLOGY OF A CARTESIAN PRODUCT OF GROUPS 

(DEFUN GR-CRTS-PRDC-FINTR (group1 group2)
  (declare
   (type group group1 group2))
  (flet ((rslt (dmns bar)
               (with-slots ((nullel1 nullel) (cmpr1 cmpr)) group1
                 (with-slots ((nullel2 nullel) (cmpr2 cmpr)) group2
                   (do ((indx 0 (1+ indx))
                        (rslt-dgop1 0)
                        (rslt-bar1 +empty-list+)
                        (rslt-dgop2 0)
                        (rslt-bar2 +empty-list+)
                        (bark 1 (ash bark +1)))
                       ((> indx (1- dmns)) 
                        (progn
                          (multiple-value-bind (dgop dgop1 dgop2) (extract-common-dgop rslt-dgop1 rslt-dgop2)
                            (declare (fixnum dgop dgop1 dgop2))
                            (absm dgop (crpr dgop1 (nreverse rslt-bar1) dgop2 (nreverse rslt-bar2))))))
                     (let ((grcrpr (pop bar)))
                       (with-grcrpr (e1 e2) grcrpr
                         (if (equal (funcall cmpr1 e1 nullel1) :equal)
                             (incf rslt-dgop1 bark)
                           (push e1 rslt-bar1))
                         (if (equal (funcall cmpr2 e2 nullel2) :equal)
                             (incf rslt-dgop2 bark)
                           (push e2 rslt-bar2)))))))))
    (the intr-mrph #'rslt)))


(DEFUN GR-CRTS-PRDC-GINTR (group1 group2)
  (declare
   (type group group1 group2))
  (flet ((rslt (dmns crpr)
               (with-slots ((nullel1 nullel)) group1
                 (with-slots ((nullel2 nullel)) group2
                   (with-crpr (dgop1 gmsm1 dgop2 gmsm2) crpr
                     (let ((bar1 (g-absm-bar nullel1 (absm dgop1 gmsm1)))
                           (bar2 (g-absm-bar nullel2 (absm dgop2 gmsm2))))
                       
                       
                       (do ((indx 0 (1+ indx))
                            (rslt-bar +empty-list+)
                            )
                           ((> indx (1- dmns)) (g-bar-absm (grcrpr nullel1 nullel2) (nreverse rslt-bar)))
                         (let ((item (grcrpr (pop bar1) (pop bar2))))
                           (push item rslt-bar)))))))))
    
    (the intr-mrph #'rslt)))


(DEFUN GR-CRTS-PRDC-RDCT (group1 group2)
  (declare
   (type group group1 group2))
  (let* ((tcc (k-g-1 (gr-crts-prdc group1 group2)))
         (bcc (crts-prdc (k-g-1 group1) (k-g-1 group2)))
         (f (build-smmr :sorc tcc :trgt bcc :degr 0 :sintr (gr-crts-prdc-fintr group1 group2)
                        :orgn `(gr-crts-prdc-rdct-f ,group1 ,group2)))
         (g (build-smmr :sorc bcc :trgt tcc :degr 0 :sintr (gr-crts-prdc-gintr group1 group2)
                        :orgn `(gr-crts-prdc-rdct-g ,group1 ,group2)))
         (h (zero-mrph tcc tcc 1)))
    (build-rdct :f f :g g :h h :orgn `(gr-crts-prdc-rdct ,group1 ,group2))))


(DEFUN GR-CRTS-PRDC-EFHM (group1 group2)
  (declare
   (type group group1 group2))
  (let* ((rdct (gr-crts-prdc-rdct group1 group2))
         (chcm1 (tcc rdct))
         (chcm2 (bcc rdct))
         (leq (build-hmeq :lrdct (trivial-rdct chcm1) :rrdct rdct :orgn `(gr-crts-prdc-leq ,group1 ,group2)))
         (req (efhm chcm2)))
    (cmps leq req)))


#|
 (cat-init)
 (setf group1 (cyclicgroup 3) group2 (cyclicgroup 2))
 (setf kg11 (k-g-1 group1) kg21 (k-g-1 group2))
 (setf g1xg2 (gr-crts-prdc group1 group2))
 (setf kg1xg2 (k-g-1 g1xg2))
 (setf kg11xkg21 (crts-prdc kg11 kg21))
 (setf rdct  (gr-crts-prdc-rdct group1 group2))
 (pre-check-rdct rdct)
 (setf *tc* (cmbn 4 1 (list (grcrpr 0 1) (grcrpr 1 0) (grcrpr 2 1) (grcrpr 2 0))))
 (setf *bc* (cmbn 4 1 (crpr 5 '(1 2) 10 '(1 1))))
 (check-rdct)
 
 (setf *tc* (cmbn 5 4 (list (grcrpr 0 1) (grcrpr 1 0) (grcrpr 2 1) (grcrpr 2 0) (grcrpr 0 1))))
 (setf *bc* (cmbn 5 1 (crpr 0 '(1 2 1 2 2) 0 '(1 1 1 1 1))))


 (check-rdct)

 (setf *bc* (cmbn 5 1 (crpr 0 '(1 2 1 2 2) 0 '(1 1 1 1 1)) -3 (crpr 5 '(1 2 2) 10 '(1 1 1)) ))
 (check-rdct)
 
 (setf efhm (gr-crts-prdc-efhm group1 group2))
 (setf (slot-value kg1xg2 'efhm) efhm)


|#     

;;;; CENTRAL EXTENSIONS


#|

(progn 
  (setf p 5 n 3)
  (setf A (cyclicGroup (expt  p (- n 2))))
  (setf C1 (cyclicGroup p))
  (setf C (gr-crts-prdc c1 c1))
  (setf KA1 (k-g-1 A))
  (setf KC1 (k-g-1 C))                    
  (setf efhm (gr-crts-prdc-efhm c1 c1))
  (setf (slot-value kc1 'efhm) efhm)
  
  
  (setf cocycle #'(lambda (crpr1 crpr2)
                    (with-grcrpr (x1 y1) crpr1
                      (with-grcrpr (x2 y2) crpr2
                        (mod (* y1 x2 (1- p) (expt p (- n 3))) (expt p (- n 2)))))))
  )
|#

(DEFUN GR-CNTR-EXTN-MULT (mult1 mult2 cocycle)
  (declare 
   (type function mult1 mult2))
  (flet ((rslt (grcrpr1 grcrpr2)
               (with-grcrpr (e11 e12) grcrpr1
                 (with-grcrpr (e21 e22) grcrpr2
                   (grcrpr (funcall mult1 (funcall mult1 e11 e21) (funcall cocycle e12 e22))
                           (funcall mult2 e12 e22))))))
    #'rslt))


(DEFUN GR-CNTR-EXTN-INV (inv1 inv2 mult1 cocycle)
  (declare 
   (type function inv1 inv2))
  (flet ((rslt (grcrpr)
               (with-grcrpr (e1 e2) grcrpr
                 (grcrpr (funcall inv1  (funcall mult1 e1 (funcall cocycle e2 (funcall inv2 e2 ))))
                         (funcall inv2 e2)))))
    #'rslt))

(DEFUN GR-CNTR-EXTN (group1 group2 cocycle)
  (declare
   (type group group1 group2))
  (with-slots ((elements1 elements) (cmpr1 cmpr) (mult1 mult) (inv1 inv) (nullel1 nullel) (orgn1 orgn)) group1
    (with-slots ((elements2 elements) (cmpr2 cmpr) (mult2 mult) (inv2 inv) (nullel2 nullel) (orgn2 orgn)) group2
      (let ((elements (GR-CRTS-PRDC-ELEMENTS elements1 elements2))
            (cmpr (gr-crts-prdc-cmpr cmpr1 cmpr2))
            (mult (gr-cntr-extn-mult mult1 mult2 cocycle))
            (inv (gr-cntr-extn-inv inv1 inv2 mult1 cocycle))
            (nullel (grcrpr nullel1 nullel2))
            (orgn `(gr-cntr-extn ,group1 ,group2 ,cocycle)))
        (build-group :elements elements :cmpr cmpr :mult mult :inv inv :nullel nullel :orgn orgn)))))

#|
(setf g (gr-cntr-extn a c cocycle))
(funcall (inv1 g)  (grcrpr 1 (grcrpr 3 2)))

(setf p 3 n 3)
(setf A (cyclicGroup (expt  p (- n 2))))
(setf C1 (cyclicGroup p))
(setf C (gr-crts-prdc c1 c1))
(setf g (gr-cntr-extn a c cocycle))

(terpri)
(mapcar #'(lambda (g1)
            (mapcar #'(lambda (g2)
                        (format t " ~D * ~D = ~D "
                          g1 g2 (funcall (mult1 g) g1 g2))
                        (terpri))
              (elements g)))
  (elements g))


(setf p 2 n 3)
(setf A (cyclicGroup (expt  p (- n 2))))
(setf C1 (cyclicGroup p))
(setf C (gr-crts-prdc c1 c1))
(setf g (gr-cntr-extn a c cocycle))


(mapcar #'(lambda (g1)
            (mapcar #'(lambda (g2)
                        (if (and (not (equal (nullel g) g1))
                                 (not (equal (nullel g) g2))
                                 (not (equal (nullel g) (funcall (mult1 g) g1 g2))))
                            (with-grcrpr (a1 c1) g1
                              (with-grcrpr (c11 c12) c1
                                (with-grcrpr (a2 c2) g2
                                  (with-grcrpr (c21 c22) c2
                                    (with-grcrpr (a3 c3) (funcall (mult1 g) g1 g2)
                                      (with-grcrpr (c31 c32) c3
                                        (progn
                                          (format t "||(~D,(~D,~D))|| <= ||(~D,(~D,~D))|| + ||(~D,(~D,~D))|| "
                                            a3 c31 c32 a1 c11 c12 a2 c21 c22)
                                          (terpri) (terpri))))))))))
              (elements g)))
  (elements g)
  
  |# 


;;; EFFECTIVE HOMOLOGY OF A CENTRAL EXTENSION
 
(DEFUN COCYCLE-FIBRATION-INTR (group1 group2 cocycle-intr)
  (declare
   (type function cocycle-intr))
  (flet ((rslt (dmns bar)
               (with-slots ((mult1 mult) (inv1 inv) (nullel1 nullel)) group1
                 (with-slots ((mult2 mult) (inv2 inv) (nullel2 nullel) (cmpr2 cmpr)) group2
                   (do ((indx 0 (1+ indx))
                        (rslt-dgop 0)
                        (rslt-bar +empty-list+)
                        (bark 1 (ash bark +1)))
                       ((> indx (- dmns 2)) (absm rslt-dgop (nreverse rslt-bar)))
                     (let ((gindx (nth indx bar))
                           (cmps2 nullel1))
                       (do ((indx2 (1+ indx) (1+ indx2)))
                           ((> indx2 (- dmns 2)) ())
                         (setf cmps2 (funcall mult1 cmps2 (nth indx2 bar))))
                       (let* ((cmps1 (funcall mult1 cmps2 (nth (1- dmns) bar)))
                              (item (funcall mult2 (funcall inv2 (funcall cocycle-intr gindx cmps1))
                                             (funcall cocycle-intr gindx cmps2))))
                         (if (equal item nullel2)
                             (incf rslt-dgop bark)
                           (push item rslt-bar)))))))))
    #'rslt))


#|                   
 (setf tau-intr (cocycle-fibration-intr c a cocycle))
 (funcall tau-intr 3 (list (grcrpr 2 1) (grcrpr 2 2) (grcrpr 1 0)))                    
 (funcall tau-intr 5 (list (grcrpr 1 1) (grcrpr 2 2) (grcrpr 1 2) (grcrpr 2 1) (grcrpr 0 2)))
|#

(DEFUN COCYCLE-FIBRATION (group1 group2 cocycle-intr)
  (declare
   (type function cocycle-intr))
  (build-smmr :sorc (k-g-1 group1) :trgt (k-g-1 group2) :degr -1 :sintr (cocycle-fibration-intr group1 group2 cocycle-intr)
              :orgn  `(central extension fibration for the cocycle ,cocycle-intr)))

#|
 (setf tau (cocycle-fibration c a cocycle))
 (setf x (fibration-total tau))
 (homology x 0 3)
|#


(DEFUN COCYCLE-FIBR-ISO1-INTR (group1 group2 cocycle-intr )
  (flet ((rslt (dmns bar)
               (with-slots ((mult1 mult) (inv1 inv) (nullel1 nullel)) group1
                 (with-slots ((mult2 mult) (inv2 inv) (nullel2 nullel) (cmpr2 cmpr)) group2
                   (do ((indx 0 (1+ indx))
                        (rslt-dgop1 0)
                        (rslt-dgop2 0)
                        (rslt-bar1 +empty-list+)
                        (rslt-bar2 +empty-list+)
                        (bark 1 (ash bark +1)))
                       ((> indx (- dmns 1)) (2ABSM-ACRPR (absm rslt-dgop1 (nreverse rslt-bar1)) (absm rslt-dgop2 (nreverse rslt-bar2))))
                     (let ((gindx (nth indx bar))                         
                           (cmps2 nullel1))
                       (with-grcrpr (a c) gindx
                         (if (equal c nullel1)
                             (incf rslt-dgop1 bark)
                           (push c rslt-bar1))
                         
                         (do ((indx2 (1+ indx) (1+ indx2)))
                             ((> indx2 (- dmns 1)) ())
                           (setf cmps2 (funcall mult1 cmps2 (grcrpr-e2 (nth indx2 bar)))))
                         (let* ((item (funcall mult2 a (funcall cocycle-intr c cmps2))))
                           
                           (if (equal item nullel2)
                               (incf rslt-dgop2 bark)
                             (push item rslt-bar2))))))))))
    #'rslt))


(DEFUN COCYCLE-FIBR-ISO1 (group1 group2 cocycle-intr)
  (declare
   (type function cocycle-intr))
  (build-smmr :sorc (k-g-1 (gr-cntr-extn group2 group1 cocycle-intr)) 
              :trgt (fibration-total (cocycle-fibration group1 group2 cocycle-intr))
              :degr 0 :sintr (cocycle-fibr-iso1-intr group1 group2 cocycle-intr)
              :orgn  `(cocycle-fibr-iso1 for the cocycle ,cocycle-intr)))


(DEFUN COCYCLE-FIBR-ISO2-INTR (group1 group2 cocycle-intr )
  (flet ((rslt (dmns crpr)
               (with-slots ((mult1 mult) (inv1 inv) (nullel1 nullel)) group1
                 (with-slots ((mult2 mult) (inv2 inv) (nullel2 nullel) (cmpr2 cmpr)) group2
                   (with-crpr (dgop1 gbar1 dgop2 gbar2) crpr
                     (let ((cbar (g-absm-bar nullel1 (absm dgop1 gbar1)))
                           (abar (g-absm-bar nullel2 (absm dgop2 gbar2))))
                       (do ((indx 0 (1+ indx))
                            (rslt-dgop 0)
                            (rslt-bar +empty-list+)
                            (bark 1 (ash bark +1)))
                           ((> indx (- dmns 1)) (absm rslt-dgop (nreverse rslt-bar)))
                         (let ((cindx (nth indx cbar))
                               (aindx (nth indx abar))
                               (cmps2 nullel1))                                              
                           (do ((indx2 (1+ indx) (1+ indx2)))
                               ((> indx2 (- dmns 1)) ())
                             (setf cmps2 (funcall mult1 cmps2 (nth indx2 cbar))))
                           (let* ((aitem (funcall mult2 aindx (funcall inv2 (funcall cocycle-intr cindx cmps2))))
                                  (item (grcrpr aitem cindx )))                                    
                             (if (equal item (grcrpr nullel2 nullel1))
                                 (incf rslt-dgop bark)
                               (push item rslt-bar)))))))))))
    #'rslt))


(DEFUN COCYCLE-FIBR-ISO2 (group1 group2 cocycle)
  (declare
   (type function cocycle))
  (build-smmr :sorc (fibration-total (cocycle-fibration group1 group2 cocycle)) 
              :trgt  (k-g-1 (gr-cntr-extn group2 group1 cocycle))
              :degr 0 :sintr (cocycle-fibr-iso2-intr group1 group2 cocycle)
              :orgn  `(cocycle-fibr-iso2 for the cocycle ,cocycle)))


(DEFUN CENTRAL-EXTENSION-EFHM (group1 group2 cocycle)
  (let* ((x (fibration-total (cocycle-fibration group2 group1 cocycle)))
         (iso1 (cocycle-fibr-iso1 group2 group1 cocycle))
         (iso2 (cocycle-fibr-iso2 group2 group1 cocycle))
         (x-efhm (efhm x))
         (rrdct (rrdct x-efhm))
         (lf (lf x-efhm))
         (lg (lg x-efhm))
         (lh (lh x-efhm))
         (lrdct (build-rdct :f (cmps iso2 lf) :g (cmps lg iso1) :h lh :orgn  `(central-extension-lrdct for ,group1 ,group2 ,cocycle))))
    (build-hmeq :lrdct lrdct :rrdct rrdct :orgn `(central-extension-efhm for ,group1 ,group2 ,cocycle))))



