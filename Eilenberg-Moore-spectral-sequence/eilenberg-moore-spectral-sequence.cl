;;;  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM
;;;  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM
;;;  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM  FIBER-SPACE-EFHM

(IN-PACKAGE #:cat)

(PROVIDE "fiber-space-efhm")

(DEFUN FIBER-HAT-U-U (fibration
                      &aux (space (sorc fibration))
                      (fiber (trgt fibration))
                      (cobar (cobar space))
                      (spac-tnsr-fiber (tnsr-prdc space fiber)))
  (declare
   (type fibration fibration)
   (type simplicial-set space)
   (type chain-complex cobar spac-tnsr-fiber)
   (type simplicial-group fiber))
  (the chain-complex
    (tnsr-prdc cobar spac-tnsr-fiber)))

(DEFUN FIBER-HAT-LEFT-PERTURBATION-INTR (fibration &aux (space (sorc fibration)))
  (declare 
   (type fibration fibration)
   (type simplicial-set space))
  (let ((cprd (dgnl space)))
    (flet ((rslt (degr tnpr)
                 (declare
                  (fixnum degr)
                  (type tnpr tnpr))
                 (with-tnpr (degr1 allp1 nil tnpr2) tnpr
                   (with-tnpr (degr21 gmsm21 degr22 loop22) tnpr2
                     (let ((cprd (cmbn-list (gnrt-? cprd degr21 gmsm21))))
                       (declare (list cprd))
                       (setf cprd (rest cprd))    ;;; because \bar{A}
                       (make-cmbn :degr (1- degr)
                                  :list
                                  (mapcar
                                      #'(lambda (term)
                                          (declare (type term term))
                                          (with-term (cffc tnpr21) term
                                            (with-tnpr (degr211 gmsm211 degr212 gmsm212) tnpr21
                                              (decf degr211)
                                              (term cffc
                                                    (tnpr
                                                     (+ degr1 degr211)
                                                     (make-allp
                                                      :list (append
                                                             (allp-list allp1)
                                                             (list (cbgn degr211 gmsm211))))
                                                     (+ degr212 degr22)
                                                     (tnpr degr212 gmsm212 degr22 loop22))))))
                                    cprd)))))))
      (the intr-mrph #'rslt))))

(DEFUN FIBER-HAT-LEFT-PERTURBATION (fibration)
  (declare (type fibration fibration))
  (the morphism
    (let ((fiber-hat-u-u (fiber-hat-u-u fibration)))
      (declare (type chain-complex fiber-hat-u-u))
      (build-mrph
       :sorc fiber-hat-u-u
       :trgt fiber-hat-u-u
       :degr -1
       :intr (fiber-hat-left-perturbation-intr fibration)
       :strt :gnrt
       :orgn `(fiber-hat-left-perturbation ,fibration)))))

(DEFUN FIBER-HAT-T-U (fibration
                      &aux 
                      (fiber-hat-u-u (fiber-hat-u-u fibration))
                      (fiber-hat-left-perturbation
                       (fiber-hat-left-perturbation fibration)))
  (declare
   (type fibration fibration)
   (type chain-complex fiber-hat-u-u)
   (type morphism fiber-hat-left-perturbation))
  (the chain-complex
    (progn
      (setf (slot-value fiber-hat-left-perturbation 'sorc) fiber-hat-u-u
        (slot-value fiber-hat-left-perturbation 'trgt) fiber-hat-u-u)
      ;; because maybe these slots have been modified when constructing
      ;;   the fiber-left-hmeq-right-reduction
      (add fiber-hat-u-u fiber-hat-left-perturbation))))


(DEFUN FIBER-SZCZARBA (fibration
                       &aux (space (sorc fibration))
                       (fiber (trgt fibration))
                       (twisted-crts-prdc (fibration-total fibration))
                       (ez (ez space fiber))
                       (fiber-dtau-d (fibration-dtau-d fibration)))
  (declare
   (type fibration fibration)
   (type simplicial-set space  twisted-crts-prdc)
   (type simplicial-group fiber)
   (type reduction ez)
   (type morphism fiber-dtau-d))
  (the (values reduction morphism)
    (multiple-value-bind (szczarba bottom-perturbation)
        (add ez fiber-dtau-d)
      (with-slots (tcc f g h) szczarba
        (setf tcc twisted-crts-prdc
          (slot-value f 'sorc) twisted-crts-prdc
          (slot-value g 'trgt) twisted-crts-prdc
          (slot-value h 'sorc) twisted-crts-prdc
          (slot-value h 'trgt) twisted-crts-prdc))
      (values szczarba bottom-perturbation))))


(DEFUN FIBER-TWISTED-TNSR-PRDC (fibration
                                &aux (szczarba (fiber-szczarba fibration)))
  (declare
   (type fibration fibration)
   (type reduction szczarba))
  (the chain-complex
    (bcc szczarba)))

(DEFUN FIBER-HAT-RIGHT-PERTURBATION (fibration
                                     &aux (space (sorc fibration))
                                     (cobar (cobar space))
                                     (fiber-hat-t-u (fiber-hat-t-u fibration)))
  (declare
   (type fibration fibration)
   (type simplicial-set space)
   (type chain-complex cobar fiber-hat-t-u))
  (the morphism
    (multiple-value-bind (szczarba bottom-perturbation)
        (fiber-szczarba space)
      (declare
       (ignore szczarba)
       (type morphism bottom-perturbation))
      (let ((rslt (tnsr-prdc (idnt-mrph cobar) bottom-perturbation)))
        (declare (type morphism rslt))
        (setf (slot-value rslt 'sorc) fiber-hat-t-u
          (slot-value rslt 'trgt) fiber-hat-t-u)
        rslt))))

(DEFUN FIBER-HAT-U-T (fibration
                      &aux (space (sorc fibration))
                      (cobar (cobar space))
                      (twisted-tnsr-prdc (fiber-twisted-tnsr-prdc fibration))
                      (fiber-hat-u-u (fiber-hat-u-u fibration)))
  (declare
   (type fibration fibration)
   (type simplicial-set space)
   (type chain-complex cobar twisted-tnsr-prdc fiber-hat-u-u))
  (the chain-complex
    (let ((rslt (tnsr-prdc cobar twisted-tnsr-prdc)))
      (declare (type chain-complex rslt))
      (setf (slot-value rslt 'grmd) fiber-hat-u-u)
      rslt)))


(DEFUN FIBER-LEFT-HMEQ-HAT (fibration
                            &aux 
                            (fiber-hat-u-t (fiber-hat-u-t fibration))
                            (fiber-hat-left-perturbation
                             (fiber-hat-left-perturbation fibration)))
  (declare
   (type fibration fibration)
   (type chain-complex fiber-hat-u-t)
   (type morphism fiber-hat-left-perturbation))
  (the chain-complex
    (add fiber-hat-u-t fiber-hat-left-perturbation)))




(DEFUN FIBER-PRE-LEFT-HMEQ-LEFT-REDUCTION-F (fibration
                                             &aux 
                                             (fiber-hat-t-u (fiber-hat-t-u fibration))
                                             (fiber (trgt fibration)))
  (declare
   (type fibration fibration)
   (type simplicial-group  fiber)
   (type chain-complex fiber-hat-t-u))      
  (the morphism
    (build-mrph  
     :sorc fiber-hat-t-u
     :trgt fiber
     :degr 0
     :intr #'fiber-pre-left-hmeq-left-reduction-intr-f
     :strt :cmbn
     :orgn `(fiber-pre-left-hmeq-left-reduction-f ,fibration))))

(DEFUN FIBER-LEFT-HMEQ-LEFT-REDUCTION-G-INTR (bbspn fbspn)
  (declare (type gmsm bbspn fbspn))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (the cmbn
                 (with-cmbn (degr list) cmbn
                   (make-cmbn :degr degr
                              :list (mapcar
                                        #'(lambda (term)
                                            (declare (type term term))
                                            (with-term (cffc loop) term
                                              (term cffc
                                                    (tnpr
                                                     0 fbspn
                                                     degr (tnpr 0 bbspn degr loop)))))
                                      list))))))
    (the intr-mrph #'rslt)))


(DEFUN FIBER-LEFT-HMEQ-LEFT-REDUCTION-G (fibration
                                         &aux (space (sorc fibration))
                                         (fiber (trgt fibration))
                                         (bbspn (bspn space))
                                         (fbspn (bspn fiber))
                                         (fiber-hat-t-u (fiber-hat-t-u fibration))
                                         )
  (declare
   (type fibration fibration)
   (type simplicial-set space fiber)
   (type gmsm bbspn fbspn)
   (type chain-complex fiber-hat-t-u))      
  (the morphism
    (build-mrph  
     :sorc fiber
     :trgt fiber-hat-t-u
     :degr 0
     :intr (fiber-left-hmeq-left-reduction-g-intr bbspn fbspn)
     :strt :cmbn
     :orgn `(fiber-left-hmeq-left-reduction-g ,fibration))))


(DEFUN FIBER-PRE-LEFT-HMEQ-LEFT-REDUCTION-H-INTR (fibration
                                                  &aux 
                                                  (fiber-hat-t-u (fiber-hat-t-u fibration))
                                                  (cmpr (cmpr fiber-hat-t-u)))
  (declare
   (type chain-complex fiber-hat-t-u)
   (type cmprf cmpr))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (the cmbn
                 (with-cmbn (degr list) cmbn
                   (let ((rslt (zero-cmbn (1+ degr))))
                     (declare (type cmbn rslt))
                     (dolist (term list)
                       (declare (type term term))
                       (with-term (cffc tnpr) term
                         (with-tnpr (degr1 allp1 degr2 tnpr2) tnpr
                           (unless (zerop degr1)
                             (with-tnpr (degr21 nil degr22 fgnrt22) tnpr2
                               (when (zerop degr21)
                                 (setf allp1 (allp-list allp1)) ;;;+++
                                 (let ((last-cbgn (car (last allp1))))
                                   (declare (type cbgn last-cbgn))
                                   (with-cbgn (degrl gmsml) last-cbgn
                                     (dstr-add-term-to-cmbn cmpr
                                                            cffc (tnpr (- degr1 degrl)
                                                                       (make-allp :list (butlast allp1))
                                                                       (+ degr2 degrl 1)
                                                                       (tnpr (1+ degrl) gmsml
                                                                             degr22 fgnrt22))
                                                            rslt)))))))))
                     rslt)))))
    (the intr-mrph #'rslt)))


(DEFUN FIBER-PRE-LEFT-HMEQ-LEFT-REDUCTION-H (fibration
                                             &aux (fiber-hat-t-u (fiber-hat-t-u fibration)))
  (declare
   (type fibration fibration)
   (type chain-complex fiber-hat-t-u))      
  (the morphism
    (build-mrph
     :sorc fiber-hat-t-u
     :trgt fiber-hat-t-u
     :degr +1
     :intr (fiber-pre-left-hmeq-left-reduction-h-intr fibration)
     :strt :cmbn
     :orgn `(fiber-pre-left-hmeq-left-reduction-h ,fibration))))


(DEFUN FIBER-PRE-LEFT-HMEQ-LEFT-REDUCTION (fibration)
  (declare (type fibration fibration))
  (the reduction
    (build-rdct
     :f (fiber-pre-left-hmeq-left-reduction-f fibration)
     :g (fiber-left-hmeq-left-reduction-g fibration)
     :h (fiber-pre-left-hmeq-left-reduction-h fibration)
     :orgn `(fiber-pre-left-hmeq-left-reduction ,fibration))))


(DEFUN FIBER-LEFT-HMEQ-LEFT-REDUCTION (fibration
                                       &aux (space (sorc fibration)) 
                                       (fiber-pre-left-hmeq-left-reduction
                                        (fiber-pre-left-hmeq-left-reduction space))
                                       (fiber-hat-right-perturbation
                                        (fiber-hat-right-perturbation space)))			            
  (declare
   (type simplicial-set space)
   (type reduction fiber-pre-left-hmeq-left-reduction)
   (type morphism fiber-hat-right-perturbation))
  (the reduction
    (progn
      (dstr-change-sorc-trgt fiber-hat-right-perturbation
                             :new-sorc (tcc fiber-pre-left-hmeq-left-reduction)
                             :new-trgt (tcc fiber-pre-left-hmeq-left-reduction))
      (let ((rslt (special-bpl fiber-pre-left-hmeq-left-reduction
                               fiber-hat-right-perturbation)))
        (declare (type reduction rslt))
        (with-slots (tcc f g h) rslt
          (setf tcc (fiber-left-hmeq-hat fibration)
            (slot-value f 'sorc) tcc
            (slot-value g 'trgt) tcc
            (slot-value h 'sorc) tcc
            (slot-value h 'trgt) tcc)
          rslt)))))



(DEFUN EILENBERG-MOORE-BICOMPLEX (fibration)
  (declare (type fibration fibration))
  (let ((rslt (fiber-left-hmeq-hat fibration)))
    (setf (slot-value rslt 'orgn) `(EILENBERG-MOORE-BICOMPLEX ,fibration))
    rslt))


(DEFUN LOOP-FBR (space)
  (declare (type simplicial-set space))
  (build-smmr
   :sorc space
   :trgt (loop-space space)
   :degr -1
   :sintr #'(lambda (dmns gmsm)
              (declare (ignore dmns))
              (absm 0 (loop3 0 gmsm 1)))
   :orgn `(total-fibration ,space)))


(DEFMETHOD SEARCH-EFHM (chcm (orgn (eql 'EILENBERG-MOORE-BICOMPLEX)))
  (declare (type chain-complex chcm))
  (let* ((fibration (second (orgn chcm)))
         (base (sorc fibration))
         (fiber (trgt fibration))
         (fiber-hat-left-perturbation
          (fiber-hat-left-perturbation fibration))
         (cobar (cobar base))
         (ttp (fiber-twisted-tnsr-prdc fibration))
         (ez (ez base fiber))
         (fiber-dtau-d (fibration-dtau-d fibration)))
    (declare
     (type fibration fibration)
     (type simplicial-set base ez)
     (type simplicial-group fiber)
     (type morphism fiber-hat-left-perturbation fiber-dtau-d)
     (type chain-complex cobar ttp))
    (setf (slot-value cobar 'efhm) (cobar (efhm base)))
    (multiple-value-bind (szczarba bottom-perturbation)
        (add ez fiber-dtau-d)
      (declare (ignore szczarba)
               (type morphism bottom-perturbation))
      (setf (slot-value ttp 'efhm) (add (tnsr-prdc (efhm base) (efhm fiber)) bottom-perturbation)))
     (the homotopy-equivalence
      (add 
       (tnsr-prdc (efhm cobar) (efhm ttp))
       fiber-hat-left-perturbation))))



(DEFVAR fib-cobar-flin)
(setf fib-cobar-flin 
  #'(lambda (degr tnpr)
      (let ((allp (gnrt1 tnpr)))      
        (declare (ignore degr)
                 (type allp allp))
        (with-allp (l) allp
          (the fixnum
            (- (length l)))))))

(defun allp-bidegree (allp)
  (declare (type allp allp))
  (the list
    (let ((p (1- (length allp)))
          (q 0))
      (declare (type fixnum p q))
      (setf q (apply #'+ (mapcar #'car (rest allp))))
      (incf q p)
      (return-from allp-bidegree (list (- p) q)))))

(DEFVAR cobar-flin)
(setf cobar-flin 
  #'(lambda (degr allp)
      (declare (ignore degr)
               (type allp allp))
      (with-allp (l) allp
        (the fixnum
          (- (length l))))))


(setf cobar-flin 
  #'(lambda (degr allp)
      (declare (type fixnum degr)
               (type allp allp))
      (let* ((allp-bidgr (allp-bidegree allp))
             (degr1 (first allp-bidgr))
             (degr2 (second allp-bidgr)))
        (if (=  degr (+ degr1 degr2))
            degr1
          0))))





;; Function that constructs the Eilenberg-Moore spectral sequence associated with a fibration
;; Returns an object of the class Spectral-sequence
(DEFUN FIBRATION-EILENBERG-MOORE-SPECTRAL-SEQUENCE (f)
  (declare (type fibration f))
  (let* ((k (eilenberg-moore-bicomplex f))
         (ecc (rbcc (efhm k))))
    (declare 
     (type chain-complex k)
     (type chain-complex ecc))
    (progn
      (change-chcm-to-flcc ecc fib-cobar-flin)
      (the spectral-sequence
        (build-ss ecc `(Fibration-Eilenberg-Moore-Spectral-Sequence ,f))))))


#|
(cat-init)
(progn
  (setf s3 (sphere 3))
  (setf k3 (chml-clss s3 3))
  (setf F3 (z-whitehead s3 k3))
  (setf X (fibration-total f3)))


(setf ss1 (FIBRATION-EILENBERG-MOORE-SPECTRAL-SEQUENCE f3))


(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss1 2 (- p) q))))


(setf bc (fltrcm ss1))
(homology bc 4)


(print-spsq-group ss1 2 0 4)
(print-spsq-group ss1 10 0 4)
(print-spsq-group ss1 2 -2 6)
(print-spsq-group ss1 3 -2 6)
(print-spsq-group ss1 4 -2 6)





(setf fibration f3)
(setf base (sorc fibration))
(setf fiber (trgt fibration))
(setf fiber-hat-left-perturbation
  (fiber-hat-left-perturbation fibration))
(setf cobar (cobar base))
(setf ttp (fiber-twisted-tnsr-prdc fibration))
(setf ez (ez base fiber))
(setf fiber-dtau-d (fibration-dtau-d fibration)))

(setf (slot-value cobar 'efhm) (cobar (efhm base)))
(multiple-value-bind (szczarba bottom-perturbation)
    (add ez fiber-dtau-d)
  
  (setf (slot-value ttp 'efhm) (add (tnsr-prdc (efhm base) (efhm fiber)) bottom-perturbation)))




|#



#|
(cat-init)
(progn
  (setf kz22 (k-z2 2))
  (setf k2 (chml-clss kz22 2))
  (setf F2 (z2-whitehead kz22 k2))
  (setf X (fibration-total f2)))


(setf ss2 (FIBRATION-EILENBERG-MOORE-SPECTRAL-SEQUENCE f2))


(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss2 2 (- p) q))))

(print-spsq-group ss1 2 0 4)
(print-spsq-group ss1 5 0 4)
(print-spsq-group ss1 2 -2 6)
(print-spsq-group ss1 3 -2 6)
(print-spsq-group ss1 4 -2 6)
|#



#|
(cat-init)
(progn
  (setf kz43 (k-zp 4 3 ))
  (setf k3 (chml-clss kz43 3))
  (setf F3 (zp-whitehead 4 kz43 k3))
  (setf X (fibration-total f3)))


(setf ss3 (FIBRATION-EILENBERG-MOORE-SPECTRAL-SEQUENCE f3))

(setf bc (fltrcm ss3))


(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss3 2 (- p) q))))

(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss3 4 (- p) q))))

|#



#|
(cat-init)
(progn
  (setf s2 (sphere 2))
  (setf k2 (chml-clss s2 2))
  (setf F2 (z-whitehead s2 k2))
  (setf X (fibration-total f2)))


(setf ss2 (FIBRATION-EILENBERG-MOORE-SPECTRAL-SEQUENCE f2))


(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss2 2 (- p) q))))


(dotimes (n 6)
  (dotimes (p (1+ n))
    (let ((q (+ n p)))
      (print-spsq-group ss2 3 (- p) q))))

(print-spsq-group ss1 2 0 4)
(print-spsq-group ss1 5 0 4)
(print-spsq-group ss1 2 -2 6)
(print-spsq-group ss1 3 -2 6)
(print-spsq-group ss1 4 -2 6)
|#

