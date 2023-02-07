(in-package :cat)

(DEFUN UNIVERSAL-COVER-ABSM-TWOP (smst group twop-edges)
  (declare
   (type simplicial-set smst)
   (type group group)
   (type function twop-edges)) 
  (labels ((twop-intr (dmns absm)
                      (declare (fixnum dmns)
                               (type absm absm))
                      (with-absm (dgop gmsm) absm
                        (declare
                         (fixnum dgop)
                         (type gnrt gmsm))
                        (if (eql 1 dmns)
                            (absm 0 (funcall twop-edges gmsm))
                          (if (> dgop 0)
                              (let* ((dgop-list (dgop-int-ext dgop))
                                     (dgop1 (first dgop-list)))
                                (declare (type list dgop-list)
                                         (fixnum dgop1))
                                (if (= dgop1 (1- dmns)) (absm (mask (- dmns 1)) 0)
                                  (let ((twop-1 (twop-intr (1- dmns) (absm (dgop-ext-int (cdr dgop-list)) gmsm)))
                                        )
                                    (declare (type absm twop-1))
                                    (absm (mask (- dmns 1)) (gmsm twop-1)))))                      
                            
                            (let* ((face-rslt (face smst 0 dmns absm))
                                   (twop (twop-intr (1- dmns) face-rslt)))
                              (declare (type absm face-rslt)
                                       (type absm twop))
                              (absm (mask (- dmns 1)) (gmsm twop))
                              ))))))
    #'twop-intr))


(DEFUN UNIVERSAL-COVER (smst group twop-edges)
  (declare
   (type simplicial-set smst)
   (type group group)
   (type function twop-edges))
  
  (let* ((twop-absm-intr (universal-cover-absm-twop smst group twop-edges))
         (twop (build-smmr
                :sorc smst
                :trgt (k-g group 0)
                :degr -1
                :sintr #'(lambda (dmns gmsm)
                           (funcall twop-absm-intr dmns (absm 0 gmsm)))
                :strt :gnrt
                :orgn `(Universal-cover ,smst))))
    (declare (type simplicial-mrph twop))
    (fibration-total twop)                            
    ))


;;; EXAMPLES
 
#|
(setf proj-plane (build-finite-ss
                  '(v
                    1 a (v v)
                    2 t (a v a))))


(setf proj-plane-twop-edges
  #'(lambda (edge)
      (if (eq edge 'a) 1)))

(setf proj-plane-universal-cover (universal-cover proj-plane (cyclicgroup 2) proj-plane-twop-edges))

(homology proj-plane-universal-cover 0)
(homology proj-plane-universal-cover 1)
(homology proj-plane-universal-cover 2)


(setf space-with-fundamental-group-z3 (build-finite-ss
                                       '(v
                                         1 a (v v) e (v v)
                                         2 t1 (a e a) t2 (a v e))))


(setf space-with-fundamental-group-z3-twop-edges
  #'(lambda (edge)
      (if (eq edge 'a) 1
        0)))


(setf space-with-fundamental-group-z3-universal-cover 
  (universal-cover space-with-fundamental-group-z3 (cyclicgroup 3) space-with-fundamental-group-z3-twop-edges))


(homology space-with-fundamental-group-z3-universal-cover 0)
(homology space-with-fundamental-group-z3-universal-cover 1)
(homology space-with-fundamental-group-z3-universal-cover 2)


(setf semiline (build-smst
                          :cmpr #'f-cmpr
                          :basis :locally-effective
                          :bspn 0 
                          :face #'(lambda (indx dmns gmsm)
                                    (if (eql 1 dmns)
                                        (if (eql 0 indx) (absm 0 (second gmsm))
                                          (absm 0 (first gmsm)))
                                      ;;(absm (mask (1- dmns)) 0)
                                      ))
                          
                          :orgn `(semiline)))

(setf semiline-f
  (build-mrph :sorc semiline :trgt (z-chcm) :degr 0
              :intr #'(lambda (dmns gmsm)
                        (if (eql 0 dmns) (cmbn 0 1 :Z-GNRT)
                          (cmbn dmns)))
              :strt :gnrt
              :orgn '(semiline-f)))


(setf semiline-g
  (build-mrph :sorc (z-chcm) :trgt semiline :degr 0
              :intr #'(lambda (dmns gmsm)
                        (if (eql 0 dmns) (cmbn 0 1 0)
                          (cmbn dmns)))
              :strt :gnrt
              :orgn '(semiline-g)))

(setf semiline-h
  (build-mrph :sorc semiline :trgt semiline :degr 1
              :intr #'(lambda (dmns gmsm)
                        (if (eql 0 dmns) 
                            (if (eql gmsm 0) (cmbn 1)
                              (let ((l (mapcar #'(lambda (i)
                                                   (list 1 (list i (1+ i))))
                                         (<a-b< 0 gmsm))))
                                (make-cmbn :degr 1 :list l))) 
                          (cmbn (1+ dmns))))
              :strt :gnrt
              :orgn '(semiline-h2)))


(setf semiline-rdct
  (build-rdct :f semiline-f
              :g semiline-g
              :h semiline-h
              :orgn '(semiline-rdct)))

(setf semiline-efhm
  (build-hmeq :lrdct (trivial-rdct semiline)
              :rrdct semiline-rdct))

(setf (slot-value semiline 'efhm) semiline-efhm)

;;(homology semiline 0)
;;(homology semiline 1)
;;(homology semiline 2)
;;(homology semiline 3)


(setf X (crts-prdc proj-plane semiline))
(setf X3 (crts-prdc space-with-fundamental-group-z3 semiline))


;;(homology X 0)
;;(homology X 1)
;;(homology X 2)
;;(homology X 3)


;;(homology X3 0)
;;(homology X3 1)
;;(homology X3 2)
;;(homology X3 3)



(setf X-twop-edges
  #'(lambda (edge)
      (with-crpr (dgop1 gmsm1 dgop2 gmsm2) edge
        (if (and (eq 0 dgop1) (eq gmsm1 'a)) 1 0))))

(setf X3-twop-edges ;; In fact this is the same function
  #'(lambda (edge)
      (with-crpr (dgop1 gmsm1 dgop2 gmsm2) edge
        (if (and (eq 0 dgop1) (eq gmsm1 'a)) 1 0))))



;;;(defun proj-plane-absm-twop (dmns absm)
;;;  (with-absm (dgop gmsm) absm
;;;    (if (eql 1 dmns)
;;;        (if (eql 0 dgop) (absm 0 1) (absm 0 0))
;;;      (if (> dgop 0)
;;;          (let* ((dgop-list (dgop-int-ext dgop))
;;;                 (dgop1 (first dgop-list)))
;;;            (if (= dgop1 (1- dmns)) (absm (mask (- dmns 1)) 0)
;;;              (let ((twop-1 (proj-plane-absm-twop (1- dmns) (absm (dgop-ext-int (cdr dgop-list)) gmsm)))
;;;                     )
;;;                (absm (mask (- dmns 1)) (gmsm twop-1)))))
;;;                
;;;                            
;;;      (let* ((face-rslt (face proj-plane 0 dmns absm))
;;;             (twop (proj-plane-absm-twop (1- dmns) face-rslt)))
;;;        (absm (mask (- dmns 1)) (gmsm twop))
;;;        )))))
;;;
;;;
;;;(defun proj-plane-absm-twop3 (dmns absm)
;;;  (with-absm (dgop gmsm) absm
;;;    (if (eql 1 dmns)
;;;        (if (and (eql 0 dgop) (eql gmsm 'a)) (absm 0 1) (absm 0 0))
;;;      (if (> dgop 0)
;;;          (let* ((dgop-list (dgop-int-ext dgop))
;;;                 (dgop1 (first dgop-list)))
;;;            (if (= dgop1 (1- dmns)) (absm (mask (- dmns 1)) 0)
;;;              (let ((twop-1 (proj-plane-absm-twop3 (1- dmns) (absm (dgop-ext-int (cdr dgop-list)) gmsm)))
;;;                     )
;;;                (absm (mask (- dmns 1)) (gmsm twop-1)))))
;;;                
;;;                            
;;;      (let* ((face-rslt (face space-with-fundamental-group-z3 0 dmns absm))
;;;             (twop (proj-plane-absm-twop3 (1- dmns) face-rslt)))
;;;        (absm (mask (- dmns 1)) (gmsm twop))
;;;        )))))


;;(cat-init)

(setf X-universal-cover (universal-cover X (cyclicgroup 2) X-twop-edges))
(setf X3-universal-cover (universal-cover X3 (cyclicgroup 3) X3-twop-edges))


(homology X-universal-cover 0)
(homology X-universal-cover 1)
(homology X-universal-cover 2)
(homology X3-universal-cover 0)
(homology X3-universal-cover 1)
(homology X3-universal-cover 2)


|#
