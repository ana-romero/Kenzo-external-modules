;; GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES
;; GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES
;; GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES   GENERALIZED-SPECTRAL-SEQUENCES


;;; GENERALIZED SPECTRAL SEQUENCES FOR FILTERED CHAIN COMPLEXES (LINEAR FILTRATION)


(DEFUN CONTAINS (l e)
  (declare (list l))
  (let ((rslt nil))
    (declare (list rslt))
    (dotimes (i (length l))
      (declare (type fixnum i))
      (if (eq e (nth i l))
          (setq rslt 't)))
    rslt))


(DEFUN G-SPSQ-NUM-MTRX (fltrcm z s p degr)
  (declare (type filtered-chain-complex fltrcm)
           (fixnum z s p degr))
  (the matrix
    (let* ((fltr-p-basis (fltrd-basis fltrcm degr p))
           (p-basis-l (length fltr-p-basis))
           (z1-list (fltr-chcm-z-gnrt-list fltrcm (- p z) p (- degr p)))
           (fltr-p-s-list (fltr-basis-gnrt-list fltrcm degr s))
           (z1-mtrx (gnrt-list-to-mtrx z1-list))
           (fltr-p-s-mtrx (gnrt-list-to-mtrx fltr-p-s-list))
           (num1-mtrx (mtrx-conc z1-mtrx fltr-p-s-mtrx))
           (num1-line-n (line-number num1-mtrx)))
      (declare
       (type fixnum p-basis-l num1-line-n)
       (type list fltr-p-basis  z1-list fltr-p-s-list)
       (type matrix z1-mtrx fltr-p-s-mtrx num1-mtrx))
      (if (eql num1-line-n p-basis-l)
          num1-mtrx 
        (let* ((num-line-n p-basis-l)
               (num-column-n (column-number num1-mtrx))
               (rslt  
                #-ACLPC (make-array (list num-line-n num-column-n)
                                    :element-type 'fixnum
                                    :initial-element 0)
                #+ACLPC
                (if (or (zerop num-line-n) (zerop num-column-n))
                    (make-array (list num-line-n num-column-n)
                                :element-type 'fixnum)
                  (make-array (list num-line-n num-column-n)
                              :element-type 'fixnum
                              :initial-element 0))))
          (declare
           (type fixnum num-line-n num-column-n)
           (type matrix rslt))
          (do ((j 0 (1+ j)))
              ((>= j num-column-n))
            (declare (type fixnum j))
            (do ((i 0 (1+ i)))
                ((>= i num1-line-n))
              (declare (type fixnum i))
              (setf (aref rslt i j) (aref num1-mtrx i j))))
          rslt)))))


(DEFUN G-SPSQ-DEN-MTRX (fltrcm s p b degr)
  (declare (type filtered-chain-complex fltrcm)
           (fixnum s p b degr))
  (the matrix
    (let* ((fltr-p-basis (fltrd-basis fltrcm degr p))
           (p-basis-l (length fltr-p-basis))
           (fltr-p-s-list (fltr-basis-gnrt-list fltrcm degr s))
           (z2-list (fltr-chcm-z-gnrt-list fltrcm (- b p) b (- (1+ degr) b)))
           (z2-mtrx (gnrt-list-to-mtrx z2-list))
           (fltr-p-s-mtrx (gnrt-list-to-mtrx fltr-p-s-list))
           (dffr-mtrx (flcc-dffr-mtrx fltrcm (1+ degr) b))
           (nil-mtrx 
            #-ACLPC (make-array (list 0 0)
                                :element-type 'fixnum
                                :initial-element 0)
            #+ACLPC (make-array (list 0 0)
                                :element-type 'fixnum))
           (dffr-z2-mtrx 
            (if (or (= 0 (array-total-size z2-mtrx)) (= 0 (array-total-size dffr-mtrx)))
                nil-mtrx
              (mtrx-prdc dffr-mtrx z2-mtrx)))
           (bnd-z2-line-n p-basis-l)
           (bnd-z2-column-n (column-number dffr-z2-mtrx))
           (bnd-z2-mtrx
            (if (or (= 0 bnd-z2-line-n) (= 0 bnd-z2-column-n))
                nil-mtrx
              (submatrix dffr-z2-mtrx 0 (1- bnd-z2-line-n) 0 (1- bnd-z2-column-n))))
           (den1-mtrx (mtrx-conc bnd-z2-mtrx fltr-p-s-mtrx))
           (den-column-n (column-number den1-mtrx))
           (den-line-n p-basis-l)
           (rslt  
            #-ACLPC (make-array (list den-line-n den-column-n)
                                :element-type 'fixnum
                                :initial-element 0)
            #+ACLPC
            (if (or (zerop den-line-n) (zerop den-column-n))
                (make-array (list den-line-n den-column-n)
                            :element-type 'fixnum)
              (make-array (list den-line-n den-column-n)
                          :element-type 'fixnum
                          :initial-element 0))))
      (declare
       (list fltr-p-basis fltr-p-s-list z2-list)
       (type fixnum p-basis-l bnd-z2-line-n bnd-z2-column-n den-column-n den-line-n)
       (type matrix z2-mtrx fltr-p-s-mtrx dffr-mtrx nil-mtrx 
             dffr-z2-mtrx bnd-z2-mtrx den1-mtrx rslt))
      (do ((j 0 (1+ j)))
          ((>= j den-column-n))
        (declare (type fixnum j))
        (do ((i 0 (1+ i)))
            ((>= i (min den-line-n (line-number den1-mtrx))))
          (declare (type fixnum i))
          (setf (aref rslt i j) (aref den1-mtrx i j))))
      rslt)))


(DEFUN G-SPSQ-BASIS-DVS (fltrcm z s p b degr)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum z s p b degr))
  (the list
    (let* ((fltr-p-basis (fltrd-basis fltrcm degr p))
           (num-mtrx (g-spsq-num-mtrx fltrcm z s p degr))
           (den-mtrx (g-spsq-den-mtrx fltrcm s p b degr)) 
           (cmbn-list nil)
           (divs nil))
      (declare
       (list fltr-p-basis cmbn-list divs)
       (type matrix num-mtrx den-mtrx))
      (progn
        (if (= 0 (array-total-size num-mtrx))
            (progn
              (setq cmbn-list nil)
              (setq divs nil))
          (if (= 0 (array-total-size den-mtrx))
              (progn
                (setq cmbn-list (gnrt-mtrx-to-list num-mtrx))
                (setq divs (make-list (length cmbn-list) :initial-element 0)))
            (let ((quotient-list (mtrx-quotient num-mtrx den-mtrx)))
              (declare (list quotient-list))
              (setq cmbn-list (first quotient-list)
                  divs (second quotient-list)))))
        (list
         (mapcar 
             #'(lambda (int-list)
                 (declare (list int-list))
                 (the cmbn
                   (let* ((cmbn (cmbn degr))
                          (cmpr (cmpr fltrcm)))
                     (declare 
                      (type cmbn cmbn)
                      (type cmprf cmpr))
                     (do ((mark1 int-list (cdr mark1))
                          (mark2 fltr-p-basis (cdr mark2)))
                         ((or (endp mark1) (endp mark2)))
                       (declare (type list mark1 mark2))
                       (if (not (= 0 (car mark1)))
                           (setq cmbn (2cmbn-add cmpr cmbn 
                                                 (cmbn degr (car mark1) (car mark2))))))
                     cmbn))) 
           cmbn-list)
         divs)))))


(DEFUN G-SPSQ-GNRTS (fltrcm z s p b degr)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum z s p b degr))
  (let* ((basis-dvs (g-spsq-basis-dvs fltrcm z s p b degr))
         (basis (first basis-dvs))
         (dvs (second basis-dvs))
         (rslt nil))
    (declare (type list basis-dvs basis dvs rslt))
    (do ((mark1 basis (cdr mark1))
         (mark2 dvs (cdr mark2)))
        ((endp mark1))
      (declare (type list mark1 mark2))
      (if (not (eql 1 (car mark2)))
          (setq rslt (append rslt (list (car mark1))))))
    rslt))


(DEFUN G-SPSQ-GROUP (fltrcm z s p b degr)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum z s p b degr))  
  (let* ((basis-dvs (g-spsq-basis-dvs fltrcm z s p b degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Generalized spectral sequence S[~D,~D,~D,~D]_{~D}" z s p b degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))



;;; NEW CLASSES AND FUNCTIONS FOR GENERALIZED FILTRATIONS


;;; PARTIAL ORDERED SETS


(DEFTYPE POCMPR () '(member :less :equal :greater :undefined))


(DEFTYPE POCMPRF () 'function)
;; (function (gnrt gnrt) poset-cmpr)


(DEFMACRO POCMPR (&rest rest)
  (ecase (length rest)
    (1 `(pocmpr1 ,@rest))
    (3 `(pocmpr3 ,@rest))))


(DEFMACRO POCMPR3 (object item1 item2)
  `(funcall (pocmpr1 ,object) ,item1 ,item2))


(DEFCLASS PARTIALLY-ORDERED-SET ()
  ((pocmpr :type pocmprf :initarg :pocmpr :reader pocmpr1)
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)
   (orgn :type list :initarg :orgn :reader orgn)
   ))


(DEFMETHOD INITIALIZE-INSTANCE :after ((poset partially-ordered-set) &rest rest)
  (set (intern (format nil "K~D" (idnm poset))) poset))


(DEFVAR *poset-list*)
(SETF *poset-list* +empty-list+)
(PUSHNEW '*poset-list* *list-list*)


(DEFMETHOD PRINT-OBJECT ((poset partially-ordered-set) stream)
  (the partially-ordered-set
    (progn
      (format stream "[K~D Partially-Ordered-Set" (idnm poset))
      poset)))


(DEFUN POSET (idnm)
  (declare (type fixnum idnm))
  (the (or partially-ordered-set null)
    (find idnm *poset-list* :key #'idnm)))


(DEFUN BUILD-POSET (&key pocmpr orgn)
  (declare 
   (type pocmprf pocmpr)
   (list orgn))
  (the partially-ordered-set
    (progn
      (let ((already (find orgn *poset-list* :test #'equal :key #'orgn)))
        (declare (type (or partially-ordered-set null) already))
        (when already
          (return-from build-poset already)))
      (let ((poset (make-instance 'partially-ordered-set
                     :pocmpr pocmpr
                     :orgn orgn)))
        (declare (type partially-ordered-set poset))   
        (push poset *poset-list*)
        poset))))


(DEFUN DOWNSET-LESS= (pcmpr l1 l2)
  (declare 
   (type pocmprf pcmpr)
   (type list l1 l2))
  (mapcar #'(lambda (p1)
              (let ((less= nil))
                (mapcar #'(lambda (p2)
                            (let ((c (funcall pcmpr p1 p2)))
                              (progn
                                (if (or (eq :equal c) (eq :less c)) (setq less= 't)))))
                  l2)
                (if (not less=) (return-from downset-less= nil))  ))
    l1)
  (return-from downset-less= 't))


(DEFUN DOWNSET-CMPR (pcmpr)
  (declare (type pocmprf pcmpr))
  (flet ((rslt (l1 l2)
               (declare (type list l1 l2))
               (let ((less= (downset-less= pcmpr l1 l2))
                     (greater= (downset-less= pcmpr l2 l1)))
                 (if (and less= greater=) :equal
                   (if less= :less
                     (if greater= :greater
                       :undefined))))))
    (the pocmprf #'rslt)))


(DEFUN DOWNSETS (pos)
  (declare (type partially-ordered-set pos))
  (let ((pcmpr (pocmpr pos)))
    (declare 
     (type pocmprf pcmpr))
    (build-poset :pocmpr (downset-cmpr pcmpr)
                   :orgn   `(downset ,pos))))


(DEFUN Z2 ()
  (build-poset :pocmpr
               #'(lambda (p1 p2)
                   (declare (type list p1 p2))
                   (let ((p11 (first p1))
                         (p12 (second p1))
                         (p21 (first p2))
                         (p22 (second p2)))
                     (declare (fixnum p11 p12 p21 p22))
                     (if (= p11 p21)
                         (if (= p12 p22) :equal
                           (if (< p12 p22) :less :greater))
                       (if (< p11 p21)
                           (if (<= p12 p22) :less
                             :undefined)
                         (if (>= p12 p22) :greater
                           :undefined)))))
               :orgn '(z2)))


(DEFUN DZ2 ()
  (downsets (z2)))



;;; GENERALIZED FITERED CHAIN COMPLEXES 


(DEFTYPE GEN-CHCM-FLIN () 'function) 
;; '(function (degr gnrt) list)


(DEFCLASS GENERALIZED-FILTERED-CHAIN-COMPLEX (chain-complex)
  ((pos :type partially-ordered-set :initarg :pos :reader pos)
   (gen-flin :type gen-chcm-flin :initarg :gen-flin :reader gen-flin1)))


(DEFVAR *gflcc-list*)
(SETF *gflcc-list* +empty-list+)
(PUSHNEW '*gflcc-list* *list-list*)


(DEFMETHOD PRINT-OBJECT ((flcm generalized-filtered-chain-complex) stream)
  (the generalized-filtered-chain-complex
    (progn
      (format stream "[K~d Generalized-Filtered-Chain-Complex]" (idnm flcm))
      flcm)))


(DEFUN GFLCC (idnm)
  (declare (type fixnum idnm))
  (the (or generalized-filtered-chain-complex null)
    (find idnm *gflcc-list* :key #'idnm)))


(DEFUN BUILD-GFLCC
    (&key cmpr basis bsgn intr-dffr dffr-strt poset gen-flin orgn)
  (declare
   (type cmprf cmpr)
   (type intr-mrph intr-dffr)
   (type basis basis)
   (type gnrt bsgn)
   (type strt dffr-strt)
   (type partially-ordered-set poset)
   (type gen-Chcm-Flin gen-flin)
   (list orgn))
  (the GENERALIZED-FILTERED-CHAIN-COMPLEX
    (progn
      (let ((already (find orgn *gflcc-list* :key #'orgn :test #'equal)))
        (declare (type (or null GENERALIZED-FILTERED-CHAIN-COMPLEX) already))
        (when already
          (return-from build-gflcc already)))
      (let ((rslt (build-chcm :cmpr cmpr :basis basis :bsgn bsgn
                              :intr-dffr intr-dffr :strt dffr-strt
                              :orgn orgn)))
        (declare (type chain-complex rslt))
        (change-class rslt 'GENERALIZED-FILTERED-CHAIN-COMPLEX)
        (setf (slot-value rslt 'pos) poset)
        (setf (slot-value rslt 'gen-flin) gen-flin)
        (push rslt *gflcc-list*)
        rslt))))


(DEFMETHOD CHANGE-CHCM-TO-GFLCC ((chcm chain-complex) poset gen-flin gen-flin-orgn)
  (declare
   (type chain-complex chcm)
   (type partially-ordered-set poset)
   (type gen-chcm-flin gen-flin)
   (list gen-flin-orgn))
  (the generalized-filtered-chain-complex
    (progn     
      (let* ((orgn (append (orgn chcm) (list 'GENERALIZED-FILTERED-WITH) gen-flin-orgn))
             (already (find orgn *gflcc-list* :key #'orgn :test #'equal)))
        (declare (type (or null GENERALIZED-FILTERED-CHAIN-COMPLEX) already))
        (when already
          (return-from change-chcm-to-gflcc already))
        (let ((rslt
               (BUILD-GFLCC 
                :cmpr (cmpr chcm) :basis (basis chcm) :bsgn (bsgn chcm) :intr-dffr (intr (dffr chcm)) :dffr-strt (strt (dffr chcm))
                :poset poset :gen-flin gen-flin :orgn orgn)))
          (push rslt *gflcc-list*)
          rslt)))))


(DEFMACRO GEN-FLIN (&rest rest)
  (ecase (length rest)
    (1 `(gen-flin1 ,@rest))
    (2 `(gen-flin2 ,@rest))
    (3 `(gen-flin3 ,@rest))))


(DEFUN GEN-FLIN3 (object degr gnrt)         
  (declare
   (type filtered-chain-complex object)
   (fixnum degr)
   (type gnrt gnrt))
  (the list
    (with-slots (gen-flin) object
      (declare (type gen-chcm-flin gen-flin))
      (funcall gen-flin degr gnrt))))

;;;(DEFUN gen-FLIN2 (object cmbn)


(DEFUN GEN-FLTRD-BASIS (gfltrcm degr gen-fltr-index)
  (declare (type generalized-filtered-chain-complex gfltrcm)
           (type fixnum degr gen-fltr-index))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "Gen-Fltrd-Basis cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the list
    (let* ((pos (pos gfltrcm))
           (basis (basis gfltrcm degr))
           (rslt +empty-list+))
      (declare 
       (type list basis rslt)
       (type partially-ordered-set pos)
       )
      (mapcar #'(lambda (gnrt)
                  (let ((gflin (gen-flin gfltrcm degr gnrt)))
                    (declare (type list gflin))
                    (mapcar #'(lambda (el)
                                (declare (type gnrt el))
                                (if (and (or (eq :equal (pocmpr pos el gen-fltr-index))
                                             (eq :less (pocmpr pos el gen-fltr-index)))
                                         (not (eq (first rslt) gnrt)))
                                    (push gnrt rslt)))
                      gflin)))
        basis)
      (nreverse rslt))))


(DEFUN GEN-FLTRD-GFLIN-ELEMENTS (gfltrcm degr gflin)
  (declare (type generalized-filtered-chain-complex gfltrcm)
           (type fixnum degr)
           (type list gflin))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "Gen-Fltrd-Gflin-Elements cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the list
    (let* ((pos (pos gfltrcm))
           (basis (basis gfltrcm degr))
           (rslt +empty-list+))
      (declare 
       (type list basis rslt)
       (type partially-ordered-set pos)
       )
      (mapcar #'(lambda (gnrt)
                  (let ((gflin2 (gen-flin gfltrcm degr gnrt)))
                    (declare (type list gflin2))
                    (if (eq :equal (pocmpr (downsets pos) gflin gflin2))
                        (push gnrt rslt))))
                  basis)
      (nreverse rslt))))


(DEFUN GEN-FLIN-TO-FLIN (gflin pocmpr path)
  (declare (type gen-chcm-flin gflin)
           (type pocmprf pocmpr)
           (type list path))
  (flet ((rslt (degr gnrt)
               (declare (fixnum degr)
                        (type gnrt gnrt))
               (the fixnum
                 (let ((gf (funcall gflin degr gnrt)))
                   (declare (type list gf))
                   (dotimes (i (length path))
                     (declare (fixnum i))
                     (let ((p (nth i path)))
                       (declare (gnrt p))
                       (if (position p gf :test #'(lambda (p1 p2)
                                                    (or
                                                     (eq (funcall pocmpr p1 p2) :equal)
                                                     (eq (funcall pocmpr p1 p2) :greater))))
                           (return-from rslt i))))))))
    (the chcm-flin #'rslt)))


(DEFUN BUILD-FILTRATION-FROM-GENERALIZED-FILTRATION (gfltrcm path)
  (declare (type generalized-filtered-chain-complex gfltrcm)
           (type list path))
  (with-slots (cmpr basis bsgn dffr orgn pos gen-flin) gfltrcm
    (declare (type cmprf cmpr)
             (type list basis orgn)
             (type gnrt bsgn)
             (type morphism dffr)
             (type partially-ordered-set pos)
             (type gen-chcm-flin gen-flin))
    (build-flcc :cmpr cmpr :basis basis :bsgn :bsgn :intr-dffr (intr dffr) :dffr-strt (strt dffr)
                :flin (gen-flin-to-flin gen-flin (pocmpr pos)  path) :orgn  `(Filtered chain complex of ,gfltrcm ,path))))


;;; GENERALIZED SPECTRAL SEQUENCES 

(DEFUN GEN-FLTR-CHCM-Z-MTRX (gfltrcm z p degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z p))
  (the matrix
    (let* ((cmpr (cmpr1 gfltrcm))
           (pcmpr (pocmpr (pos gfltrcm)))
           (dffr (dffr gfltrcm))
           (sbasis (gen-fltrd-basis gfltrcm degr p))
           (tbasis1 (gen-fltrd-basis gfltrcm (1- degr) p))
           (tbasis nil))
      (declare
       (type cmprf cmpr)
       (type pocmprf pcmpr)
       (type morphism dffr)
       (type list sbasis tbasis1 tbasis))
      (progn
        (mapcar 
            #'(lambda (gnrt)
                (declare (type gnrt gnrt))
                (let ((gflin (gen-flin gfltrcm (1- degr) gnrt))
                      (b 't))
                  (declare (type list gflin))
                  (mapcar #'(lambda (el)
                              (if (or (eq :equal (funcall pcmpr el z))
                                      (eq :less (funcall pcmpr el z)))
                                  (setq b nil)))
                    gflin)
                  (if b (push gnrt tbasis))))
          tbasis1)
        (setq tbasis (nreverse tbasis))   
        (let ((srank (length sbasis))
              (trank (length tbasis)))
          (declare (type fixnum srank trank))
          (let ((rslt 
                 #-ACLPC
                 (make-array (list trank srank)
                             :element-type 'fixnum
                             :initial-element 0)
                 #+ACLPC
                 (if (or (zerop srank) (zerop trank))
                     (make-array (list trank srank)
                                 :element-type 'fixnum)
                   (make-array (list trank srank)
                               :element-type 'fixnum
                               :initial-element 0))))               
            (declare (type matrix rslt))
            (do ((j 0 (1+ j))
                 (mark sbasis (cdr mark)))
                ((endp mark))
              (declare
               (type fixnum j)
               (list mark))
              (let ((cmbn (gnrt-? dffr degr (car mark))))
                (declare (type cmbn cmbn))
                (do ((mark1 (cmbn-list cmbn) (cdr mark1)))
                    ((endp mark1))
                  (declare (list mark1))
                  (with--term (cffc gnrt) mark1
                    (declare 
                     (type fixnum cffc)
                     (type gnrt gnrt))
                    (do ((mark2 tbasis (cdr mark2))
                         (i 0 (1+ i))
                         (found nil))
                        ((or (endp mark2) found))
                      (declare 
                       (list mark2)
                       (type fixnum i found))
                      (if (eq :equal (funcall cmpr gnrt (car mark2)))
                          (progn
                            (setq found 1)
                            (setf (aref rslt i j) cffc))))))))
            rslt))))))


(DEFUN GEN-FLTR-BASIS-GNRT-LIST (gfltrcm degr fltr-index)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt fltr-index))
  (the list
    (let* ((all-basis (gen-fltrd-basis gfltrcm degr fltr-index))
           (basis-l (length all-basis)))
      (declare (type list all-basis)
               (fixnum basis-l))
            (the list
              (mapcar
                  #'(lambda (i)
                      (nconc (make-list (1- i) :initial-element 0) 
                             (list 1) (make-list (- basis-l i) :initial-element 0)))
                (<a-b> 1 basis-l))))))


(DEFUN GEN-FLTR-CHCM-Z-GNRT-LIST (gfltrcm z p degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z p))
  (the list
    (let* ((all-basis-gnrt-list (gen-fltr-basis-gnrt-list gfltrcm degr p)))
      (declare
       (list all-basis-gnrt-list))
      (if (eq :equal (pocmpr (pos gfltrcm) z p))
          all-basis-gnrt-list
        (let* ((mat (gen-FLTR-CHCM-Z-MTRX gfltrcm z p degr))
               (line-n (line-number mat))
               (column-n (column-number mat)))
          (declare 
           (type matrix mat)
           (type fixnum line-n column-n))
          (if (= 0 column-n) nil
            (if (= 0 line-n)
                all-basis-gnrt-list
              (kernel mat))))))))


(DEFUN GEN-FLTR-BASIS-2INDX-GNRT-LIST (gfltrcm degr fltr-index1 fltr-index2)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt fltr-index1 fltr-index2))
   (the list
     (let* ((basis1 (gen-fltrd-basis gfltrcm degr fltr-index1))
            (basis2 (gen-fltrd-basis gfltrcm degr fltr-index2))
            (basis-l (length basis1)))
       (declare (type list basis1 basis2)
                (type fixnum basis-l))
       (the list
         (mapcar
             #'(lambda (gnrt)
                 (let ((i (position gnrt basis1 :test #'(lambda (gnrt1 gnrt2)
                                                          (eq (cmpr gfltrcm gnrt1 gnrt2) :equal)
                                                          ))))
                   (declare (fixnum i))
                   (nconc (make-list i :initial-element 0) 
                          (list 1) (make-list (1- (- basis-l i)) :initial-element 0))))
           basis2)))))


(DEFUN GEN-SPSQ-NUM-MTRX (gfltrcm z s p degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p))
  (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         (z1-list (gen-fltr-chcm-z-gnrt-list gfltrcm z p degr))
         (fltr-s-list (gen-fltr-basis-2indx-gnrt-list gfltrcm degr p s))
         (z1-mtrx (gnrt-list-to-mtrx z1-list))
         (fltr-s-mtrx (gnrt-list-to-mtrx fltr-s-list))
         (num1-mtrx (mtrx-conc z1-mtrx fltr-s-mtrx))
         (num1-line-n (line-number num1-mtrx)))
    (declare
     (type fixnum p-basis-l num1-line-n)
     (type list fltr-p-basis z1-list fltr-s-list)
     (type matrix z1-mtrx fltr-s-mtrx num1-mtrx))
    (the matrix
      (if (eql num1-line-n p-basis-l)
          num1-mtrx 
        (let* ((num-line-n p-basis-l)
               (num-column-n (column-number num1-mtrx))
               (rslt  
                #-ACLPC (make-array (list num-line-n num-column-n)
                                    :element-type 'fixnum
                                    :initial-element 0)
                #+ACLPC
                (if (or (zerop num-line-n) (zerop num-column-n))
                    (make-array (list num-line-n num-column-n)
                                :element-type 'fixnum)
                  (make-array (list num-line-n num-column-n)
                              :element-type 'fixnum
                              :initial-element 0))))
          (declare
           (type fixnum num-line-n num-column-n)
           (type matrix rslt))
          (do ((j 0 (1+ j)))
              ((>= j num-column-n))
            (declare (type fixnum j))
            (do ((i 0 (1+ i)))
                ((>= i num1-line-n))
              (declare (type fixnum i))
              (setf (aref rslt i j) (aref num1-mtrx i j))))
          rslt)))))


(DEFUN GEN-FLCC-DFFR-MTRX (gfltrcm degr fltr-index1 fltr-index2)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (fixnum degr)
   (type gnrt fltr-index1 fltr-index2))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "gen-flcc-dffr-mtrx cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the matrix
    (let* ((cmpr (cmpr1 gfltrcm))
           (dffr (dffr gfltrcm))
           (sbasis (gen-fltrd-basis gfltrcm degr fltr-index1))
           (tbasis (gen-fltrd-basis gfltrcm (1- degr) fltr-index2))
           (srank (length sbasis))
           (trank (length tbasis)))
      (declare
       (type cmprf cmpr)
       (type morphism dffr)
       (list sbasis tbasis)
       (fixnum srank trank))
      (let ((rslt 
             #-ACLPC
             (make-array (list trank srank)
                         :element-type 'fixnum
                         :initial-element 0)
             #+ACLPC
             (if (or (zerop srank) (zerop trank))
                 (make-array (list trank srank)
                             :element-type 'fixnum)
               (make-array (list trank srank)
                           :element-type 'fixnum
                           :initial-element 0))
             ))               
        (declare (type matrix rslt))
        (do ((j 0 (1+ j))
             (mark sbasis (cdr mark)))
            ((endp mark))
          (declare
           (fixnum j)
           (list mark))
          (let ((cmbn (gnrt-? dffr degr (car mark))))
            (declare (type cmbn cmbn))
            (do ((mark1 (cmbn-list cmbn) (cdr mark1)))
                ((endp mark1))
              (declare (list mark1))
              (with--term (cffc gnrt) mark1
                (declare 
                 (fixnum cffc)
                 (type gnrt gnrt))
                (do ((mark2 tbasis (cdr mark2))
                     (i 0 (1+ i))
                     (found nil))
                    ((or (endp mark2) found))
                  (declare 
                   (list mark2)
                   (fixnum i found))
                  (if (eq :equal (funcall cmpr gnrt (car mark2)))
                      (progn
                        (setq found 1)
                        (setf (aref rslt i j) cffc))))))))
        rslt))))


(DEFUN GEN-SPSQ-DEN-MTRX (gfltrcm s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s p b))
  (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         (fltr-s-list (gen-fltr-basis-2indx-gnrt-list gfltrcm degr p s))
         (z2-list (gen-fltr-chcm-z-gnrt-list gfltrcm p b (1+ degr)))
         (z2-mtrx (gnrt-list-to-mtrx z2-list))
         (fltr-s-mtrx (gnrt-list-to-mtrx fltr-s-list))
         (dffr-mtrx (gen-flcc-dffr-mtrx gfltrcm (1+ degr) b p))
         (nil-mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (dffr-z2-mtrx 
          (if (or (= 0 (array-total-size z2-mtrx)) (= 0 (array-total-size dffr-mtrx)))
              nil-mtrx
            (mtrx-prdc dffr-mtrx z2-mtrx)))
         (bnd-z2-mtrx
          (if (or (= 0 (line-number dffr-z2-mtrx)) (= 0 (column-number dffr-z2-mtrx)))
              nil-mtrx
            dffr-z2-mtrx))
         (den1-mtrx (mtrx-conc bnd-z2-mtrx fltr-s-mtrx))
         (den-column-n (column-number den1-mtrx))
         (den-line-n p-basis-l)
         (rslt  
          #-ACLPC (make-array (list den-line-n den-column-n)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC
          (if (or (zerop den-line-n) (zerop den-column-n))
              (make-array (list den-line-n den-column-n)
                          :element-type 'fixnum)
            (make-array (list den-line-n den-column-n)
                        :element-type 'fixnum
                        :initial-element 0))))
    (declare
     (type list fltr-p-basis fltr-s-list z2-list)
     (type fixnum p-basis-l den-column-n den-line-n)
     (type matrix z2-mtrx fltr-s-mtrx dffr-mtrx nil-mtrx 
           dffr-z2-mtrx bnd-z2-mtrx den1-mtrx rslt))
    (the matrix
      (progn
      (do ((j 0 (1+ j)))
          ((>= j den-column-n))
        (declare (type fixnum j))
        (do ((i 0 (1+ i)))
            ((>= i (min den-line-n (line-number den1-mtrx))))
          (declare (type fixnum i))
          (setf (aref rslt i j) (aref den1-mtrx i j))))
      rslt))))


(DEFUN GEN-EFF-SPSQ-BASIS-DVS (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))
  (the list
    (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
           (num-mtrx (gen-spsq-num-mtrx gfltrcm z s p degr))
           (den-mtrx (gen-spsq-den-mtrx gfltrcm s p b degr)) 
           (cmbn-list nil)
           (divs nil))
      (declare
       (list fltr-p-basis cmbn-list divs)
       (type matrix num-mtrx den-mtrx))
      (progn
        (if (= 0 (array-total-size num-mtrx))
            (progn
              (setq cmbn-list nil)
              (setq divs nil))
          (if (= 0 (array-total-size den-mtrx))
              (progn
                (setq cmbn-list (gnrt-mtrx-to-list num-mtrx))
                (setq divs (make-list (length cmbn-list) :initial-element 0)))
            (let ((quotient-list (mtrx-quotient num-mtrx den-mtrx)))
              (declare (list quotient-list))
              (setq cmbn-list (first quotient-list)
                  divs (second quotient-list)))))
        (list
         (mapcar 
             #'(lambda (int-list)
                 (declare (list int-list))
                 (the cmbn
                   (let* ((cmbn (cmbn degr))
                          (cmpr (cmpr gfltrcm)))
                     (declare 
                      (type cmbn cmbn)
                      (type cmprf cmpr))
                     (do ((mark1 int-list (cdr mark1))
                          (mark2 fltr-p-basis (cdr mark2)))
                         ((or (endp mark1) (endp mark2)))
                       (declare (type list mark1 mark2))
                       (if (not (= 0 (car mark1)))
                           (setq cmbn (2cmbn-add cmpr cmbn 
                                                 (cmbn degr (car mark1) (car mark2))))))
                     cmbn))) 
           cmbn-list)
         divs)))))


(DEFUN GEN-EFF-SPSQ-GNRTS (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))
  (let* ((basis-dvs (gen-eff-spsq-basis-dvs gFltrCm z s p b degr))
         (basis (first basis-dvs))
         (dvs (second basis-dvs))
         (rslt nil))
    (declare (type list basis-dvs basis dvs rslt))
    (do ((mark1 basis (cdr mark1))
         (mark2 dvs (cdr mark2)))
        ((endp mark1))
      (declare (type list mark1 mark2))
      (if (not (eql 1 (car mark2)))
          (setq rslt (append rslt (list (car mark1))))))
    rslt))


(DEFUN GEN-EFF-SPSQ-GROUP (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))  
  (let* ((basis-dvs (gen-eff-spsq-basis-dvs gfltrcm z s p b degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Generalized spectral sequence S[~D,~D,~D,~D]_{~D}" z s p b degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))


(DEFUN GEN-EFF-SPSQ-DFFR (gfltrcm z1 s1 p1 b1 z2 s2 p2 b2 degr int-list)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z1 s1 p1 b1 z2 s2 p2 b2)
   (type list int-list))
  (the list
    (let* ((cmpr (cmpr1 gfltrcm))
           (sbasis (gen-fltrd-basis gfltrcm degr p1))
           (tbasis (gen-fltrd-basis gfltrcm (1- degr) p2))
           (num-mtrx (gen-spsq-num-mtrx gfltrcm z1 s1 p1 degr))
           (den-mtrx (gen-spsq-den-mtrx gfltrcm s1 p1 b1 degr))
           (z-gnrt-list (gen-fltr-chcm-z-gnrt-list gfltrcm  z1 p1 degr))
           (s-basis-dvs (mtrx-quotient (copy-mtrx num-mtrx) den-mtrx))
           (s-gnrt-list (first s-basis-dvs))
           (s-dvs (second s-basis-dvs))
           (t-basis-dvs (gen-eff-spsq-basis-dvs gfltrcm z2 s2 p2 b2 (1- degr)))
           (t-gnrt-list (first t-basis-dvs))
           (t-dvs (second t-basis-dvs))
           (selement (make-list (length sbasis) :initial-element 0))
           (cmbn (cmbn degr)))
      (declare 
       (type cmprf cmpr)
       (type list sbasis tbasis z-gnrt-list s-basis-dvs s-gnrt-list s-dvs
             t-basis-dvs t-gnrt-list t-dvs selement)
       (type matrix num-mtrx den-mtrx)
       (type cmbn cmbn))
      (progn
        (labels ((2list-add (list1 list2)
                            (declare (type list list1 list2))
                            (the list
                              (mapcar #'+ list1 list2)))
                 (n-list (n list)
                         (declare
                          (type list list)
                          (type fixnum n))
                         (the list
                           (mapcar #'(lambda (i)
                                       (* n i))
                             list))))
          (do ((mark1 s-gnrt-list (cdr mark1))
               (mark2 s-dvs (cdr mark2))
               (mark3 int-list))
              ((endp mark1))
            (declare (list mark1 mark2 mark3))
            (if (not (eql 1 (car mark2)))
                (progn
                  (setq selement (2list-add selement (n-list (car mark3) (car mark1))))
                  (pop mark3)))))
        (let* ((selt-coord (vctr-coordinates selement num-mtrx)))
          (declare (list selt-coord))
          (do ((mark1 z-gnrt-list (cdr mark1))
               (mark2 selt-coord (cdr mark2)))
              ((endp mark1))
            (declare (list mark1 mark2))
            (if (not (eql 0 (car mark2)))
                (setq cmbn (2cmbn-add cmpr cmbn (n-cmbn (car mark2)
                                                        (funcall 
                                                         #'(lambda (int-list)
                                                             (declare (list int-list))
                                                             (the cmbn
                                                               (let* ((cmbn2 (cmbn degr)))
                                                                 (declare 
                                                                  (type cmbn cmbn))
                                                                 (do ((mark3 int-list (cdr mark3))
                                                                      (mark4 sbasis (cdr mark4)))
                                                                     ((or (endp mark3) (endp mark4)))
                                                                   (declare (type list mark3 mark4))
                                                                   (if (not (= 0 (car mark3)))
                                                                       (setq cmbn2 (2cmbn-add cmpr cmbn2 
                                                                                              (cmbn degr (car mark3) (car mark4))))))
                                                                 cmbn2)))
                                                         (car mark1))))))))
        (let* ((dffr-cmbn (dffr gfltrcm cmbn))
               (crdnt (cmbn-coordinates dffr-cmbn t-gnrt-list tbasis cmpr))
               (rslt nil))
          (declare
           (type cmbn dffr-cmbn)
           (type list crdnt rslt))
          (do ((mark1 t-dvs (cdr mark1))
               (mark2 crdnt (cdr mark2)))
              ((endp mark1))
            (declare (list mark1 mark2))
            (if (not (eql 1 (car mark1)))
                (if (eql 0 (car mark1))
                    (setq rslt (nconc rslt (list (car mark2))))
                  (setq rslt (nconc rslt (list (mod (car mark2) (car mark1))))))))
          rslt)))))



;; USING EFFECTIVE HOMOLOGY TO DETERMINE GENERALIZED SPECTRAL SEQUENCES

(DEFUN GFLTRCM-EFCC-FLTR (gfltrcm efcc)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm efcc))
  (let* ((efhm (efhm gfltrcm)))
    (declare
     (type homotopy-equivalence efhm))
    (setf (slot-value efhm 'rbcc) efcc)))



(DEFUN GEN-SPSQ-BASIS-DVS (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))
   (let* ((hmtp-eq (efhm gfltrcm))
          (eff-chcm (rbcc hmtp-eq)))
      (declare 
        (type homotopy-equivalence hmtp-eq)
       (type chain-complex eff-chcm))
     (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
       (error "gen-spsq-basis-dvs cannot work with a locally-effective chain complex such that its effective chain complex is not a generalized filtered chain complex."))
     (progn
       (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
           (gen-eff-spsq-basis-dvs gfltrcm z s p b degr)
         (let ((spsq-list (gen-eff-spsq-basis-dvs eff-chcm z s p b degr)))
           (declare (list spsq-list))
           (if (eql gfltrcm eff-chcm) spsq-list
             (list
              (mapcar
                  #'(lambda (cmbn)
                      (lf hmtp-eq (rg hmtp-eq cmbn)))
                (first spsq-list))
              (second spsq-list))))))))
  

(DEFUN GEN-SPSQ-GNRTS (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))
   (let* ((hmtp-eq (efhm gfltrcm))
          (eff-chcm (rbcc hmtp-eq)))
      (declare 
       (type homotopy-equivalence hmtp-eq)
       (type chain-complex eff-chcm))
     (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
       (error "gen-spsq-gnrts cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
     (progn
       (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
           (gen-eff-spsq-gnrts gfltrcm z s p b degr)
         (let ((gnrts (gen-eff-spsq-gnrts eff-chcm z s p b degr)))
           (declare (list gnrts))
           (if (eql gfltrcm eff-chcm) gnrts
             (mapcar
                  #'(lambda (cmbn)
                      (lf hmtp-eq (rg hmtp-eq cmbn)))
               gnrts)))))))


(DEFUN GEN-SPSQ-GROUP (gfltrcm z s p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z s p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "gen-spsq-group cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
       (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
           (gen-eff-spsq-group gfltrcm z s p b degr)
         (gen-eff-spsq-group eff-chcm z s p b degr))))


(DEFUN GEN-SPSQ-DFFR (gfltrcm z1 s1 p1 b1 z2 s2 p2 b2 degr int-list)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z1 s1 p1 b1 z2 s2 p2 b2)
   (type list int-list))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "gen-spsq-dffr cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
           (gen-eff-spsq-dffr gfltrcm z1 s1 p1 b1 z2 s2 p2 b2 degr int-list)
      (gen-eff-spsq-dffr eff-chcm z1 s1 p1 b1 z2 s2 p2 b2 degr int-list))))





;;; FUNCTIONS FOR THE GENERALIZED SERRE SPECTRAL SEQUENCE


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
    (declare
     (type fixnum p1 p2)
     (type list rslt))
    (the list
    (progn
      (mapcar #'(lambda (i)
                  (push (2-points-add p (list i (- i))) rslt))
        (nreverse (<a-b> 1 p2)))
      (push p rslt)
      (mapcar #'(lambda (j)
                  (push (2-points-add p (list (1- (- j)) j)) rslt))
        (<a-b> 1 (1- p1)))
      rslt))))


(DEFUN E2-GSPSQ-BASIS-DVS (gfltrcm p1 p2 degr)
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


(DEFUN E2-GSPSQ-GROUP (gfltrcm p1 p2 degr)
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


(DEFUN E2-GSPSQ-GNRTS (gfltrcm p1 p2 degr)
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








