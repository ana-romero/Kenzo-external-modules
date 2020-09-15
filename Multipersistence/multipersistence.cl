;; MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE 
;; MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE 
;; MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE MULTIPERSISTENCE

(IN-PACKAGE #:cat)

(PROVIDE "multipersistence")


(DEFUN GEN-FLTR-CHCM-Z-PRST-MTRX (gfltrcm p degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p))
  (the matrix
    (let* ((cmpr (cmpr1 gfltrcm))
           (dffr (dffr gfltrcm))
           (sbasis (gen-fltrd-basis gfltrcm degr p))
           (tbasis (gen-fltrd-basis gfltrcm (1- degr) p)))
      (declare
       (type cmprf cmpr)
       (type morphism dffr)
       (type list sbasis  tbasis))
      (progn
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


(DEFUN GEN-FLTR-CHCM-Z-PRST-GNRT-LIST (gfltrcm p degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p))
  (the list
    (let* ((all-basis-gnrt-list (gen-fltr-basis-gnrt-list gfltrcm degr p)))
      (declare
       (list all-basis-gnrt-list))
      (let* ((mat (gen-fltr-chcm-z-prst-mtrx gfltrcm p degr))
             (line-n (line-number mat))
             (column-n (column-number mat)))
        (declare 
         (type matrix mat)
         (type fixnum line-n column-n))
        (if (= 0 column-n) nil
          (if (= 0 line-n)
              all-basis-gnrt-list
            (kernel mat)))))))


(DEFUN MULTIPRST-NUM-MTRX (gfltrcm p degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p))
  (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         (z1-list (gen-fltr-chcm-z-prst-gnrt-list gfltrcm p degr))
         (num1-mtrx (gnrt-list-to-mtrx z1-list))
         (num1-line-n (line-number num1-mtrx)))
    (declare
     (type fixnum p-basis-l num1-line-n)
     (type list fltr-p-basis z1-list )
     (type matrix   num1-mtrx))
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


(DEFUN MULTIPRST-DEN-MTRX (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         (z2-list (gen-fltr-chcm-z-gnrt-list gfltrcm p b (1+ degr)))
         (z2-mtrx (gnrt-list-to-mtrx z2-list))
         
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
         (den1-mtrx
          (if (or (= 0 (line-number dffr-z2-mtrx)) (= 0 (column-number dffr-z2-mtrx)))
              nil-mtrx
            dffr-z2-mtrx))
         
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
     (type list fltr-p-basis z2-list)
     (type fixnum p-basis-l den-column-n den-line-n)
     (type matrix z2-mtrx  dffr-mtrx nil-mtrx 
           dffr-z2-mtrx  den1-mtrx rslt))
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


(DEFUN EFF-MULTIPRST-BASIS-DVS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (the list
    (let* ((fltr-p-basis (gen-fltrd-basis gfltrcm degr p))
           (num-mtrx (multiprst-num-mtrx gfltrcm p degr))
           (den-mtrx (multiprst-den-mtrx gfltrcm p b degr)) 
           (cmbn-list nil)
           (divs nil))
      (declare
       (list fltr-p-basis cmbn-list divs)
       ;;(type matrix num-mtrx den-mtrx)
       )
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


(DEFUN EFF-MULTIPRST-GNRTS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((basis-dvs (EFF-MULTIPRST-basis-dvs gFltrCm p b degr))
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


(DEFUN EFF-MULTIPRST-GROUP (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))  
  (let* ((basis-dvs (EFF-MULTIPRST-basis-dvs gfltrcm p b degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Multipersistence group H[~D,~D]_{~D}" p b degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))


;; USING EFFECTIVE HOMOLOGY TO DETERMINE THE RANK INVARIANT


(DEFUN MULTIPRST-BASIS-DVS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-basis-dvs cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-multiprst-basis-dvs gfltrcm p b degr)
        (let ((m-list (eff-multiprst-basis-dvs eff-chcm p b degr)))
          (declare (list m-list))
          (if (eql gfltrcm eff-chcm) m-list
            (list
             (mapcar
                 #'(lambda (cmbn)
                     (lf hmtp-eq (rg hmtp-eq cmbn)))
               (first m-list))
             (second m-list))))))))


(DEFUN MULTIPRST-GNRTS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-gnrts cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-multiprst-gnrts gfltrcm p b degr)
        (let ((gnrts (eff-multiprst-gnrts eff-chcm p b degr)))
          (declare (list gnrts))
          (if (eql gfltrcm eff-chcm) gnrts
            (mapcar
                #'(lambda (cmbn)
                    (lf hmtp-eq (rg hmtp-eq cmbn)))
              gnrts)))))))


(DEFUN MULTIPRST-GROUP (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-group cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
        (eff-multiprst-group gfltrcm p b degr)
      (eff-multiprst-group eff-chcm p b degr))))


;;; DEFINING A NEW INVARIANT M


(DEFUN DZM-GEN-FLTRD-BASIS (gfltrcm degr ds)
  (declare (type generalized-filtered-chain-complex gfltrcm)
           (type fixnum degr )
           (type list ds))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "dzm-gen-fltrd-basis cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the list
    (let* ((pos (pos gfltrcm))
           (basis (basis gfltrcm degr))
           (rslt +empty-list+))
      (declare 
       (type list basis rslt)
       (type partially-ordered-set pos))
      (mapcar #'(lambda (gnrt)
                  (let ((gflin (gen-flin gfltrcm degr gnrt))
                        (b 't))
                    (declare (type list gflin)
                             (type gnrt gnrt))
                    (progn
                      (mapcar #'(lambda (p)
                                  (let ((bi nil))
                                    (declare (type list p))
                                    (mapcar #'(lambda (el)
                                                (if (or  
                                                     (eq :less (pocmpr pos el (list p)))
                                                     (eq :equal (pocmpr pos el (list p)))
                                                     )
                                                    (setf bi 't)))
                                      gflin)
                                    (if (not bi) 
                                        (setf b nil))))
                        ds)
                      (if (and b (not (eq (first rslt) gnrt)))
                          (push gnrt rslt)))))
        basis)
      (nreverse rslt))))


(DEFUN DZM-GEN-FLTR-CHCM-Z-MTRX (gfltrcm z p degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt z p))
  (the matrix
    (let* ((cmpr (cmpr1 gfltrcm))
           (pcmpr (pocmpr (pos gfltrcm)))
           (dffr (dffr gfltrcm))
           (sbasis (dzm-gen-fltrd-basis gfltrcm degr p))
           (tbasis1 (dzm-gen-fltrd-basis gfltrcm (1- degr) p))
           (tbasis nil)
           )
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
                  (mapcar #'(lambda (zi)
                              (let ((bi nil))
                                (declare (type list zi))
                                (mapcar #'(lambda (el)
                                            (if (or (eq :equal (funcall pcmpr el (list zi)))
                                                    (eq :less (funcall pcmpr el (list zi))))
                                                (setf bi 't)))
                                      gflin)
                                (if (not bi) 
                                    (setf b nil))))
                    z)
                    
                    
                    (if (not b) (push gnrt tbasis))))
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


(DEFUN DZM-GEN-FLTR-CHCM-Z-GNRT-LIST (gfltrcm z p degr)
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
        (let* ((mat (dzm-gen-fltr-chcm-z-mtrx gfltrcm z p degr))
               (line-n (line-number mat))
               (column-n (column-number mat)))
          (declare 
           (type matrix mat)
           (type fixnum line-n column-n))
          (if (= 0 column-n) nil
            (if (= 0 line-n)
                all-basis-gnrt-list
              (kernel mat))))))))


(DEFUN DZM-GEN-FLCC-DFFR-MTRX (gfltrcm degr fltr-index1 fltr-index2)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (fixnum degr)
   (type gnrt fltr-index1 fltr-index2))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "dzm-gen-flcc-dffr-mtrx cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the matrix
    (let* ((cmpr (cmpr1 gfltrcm))
           (dffr (dffr gfltrcm))
           (sbasis (dzm-gen-fltrd-basis gfltrcm degr fltr-index1))
           (tbasis (dzm-gen-fltrd-basis gfltrcm (1- degr) fltr-index2))
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


(DEFUN DOWNSET-LESS (d)
  (declare (type list d))
  (let ((rslt nil))
    (declare (type list rslt))
    (dolist (k (<a-b> 0 (- (length d) 1)))
      (let* ((dk (nth k d))
             (dk1 (first dk))
             (dk2 (second dk)))
        (if (and (= k 0) (> dk1 1))
            (push (list (1- dk1) dk2) rslt))
        (if (> k 0)            
            (let* ((dk-1 (nth (- k 1) d)) 
                   (dk-11 (first dk-1)))
              (if (> (- dk1 dk-11) 1) (push (list (1- dk1) dk2) rslt))))          
        (if (< k (- (length d) 1))            
            (let* ((dk+1 (nth (+ 1 k) d)) 
                   (dk+11 (first dk+1))
                   (dk+12 (second dk+1)))
              (if (or (> (- dk2 dk+12) 1) 
                      (and (= (- dk2 dk+12) 1) (= (- dk+11 dk1) 1)))                   
                  (push (list dk1 (1- dk2)) rslt))))
        (if (= k (- (length d) 1))
            (if (> dk2 1) (push (list dk1 (1- dk2)) rslt)))) )                    
    (nreverse rslt)))


(DEFUN DOWNSET-NON-COMPARABLE (d t1 t2)
  (declare (type list d)
           (fixnum t1 t2))
  (let ((rslt nil))
    (declare (type list rslt))
    (progn
      (let* ((d1 (first d))
             (d11 (first d1)))
        (declare (type list d1)
                 (fixnum d11))
        (if (> d11 1) (push (list (1- d11) t2) rslt)))
      (dolist (k (<a-b> 1 (1- (length d))))
        (let* ((dk (nth (1- k) d))
               (dk1 (first dk))
               (dk2 (second dk))
               (dk+1 (nth k d))
               (dk+11 (first dk+1))
               (dk+12 (second dk+1))
               (dk2-1 (1- dk2))
               (dk+11-1 (1- dk+11)))
          (if (and (> dk+11-1 dk1) (> dk2-1 dk+12))
              (push (list dk+11-1 dk2-1) rslt))))
      (let* ((dl (first (last d)))
             (dl2 (second dl)))
        (if (> dl2 1)
            (push (list t1 (1- dl2)) rslt))))              
      (nreverse rslt)))


(DEFUN MULTIPRST-M-NUM-MTRX-1 (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))         
         (z2-list (dzm-gen-fltr-chcm-z-gnrt-list gfltrcm p b (1+ degr)))
         (z2-mtrx (gnrt-list-to-mtrx z2-list))         
         (dffr-mtrx (dzm-gen-flcc-dffr-mtrx gfltrcm (1+ degr) b p))
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
         (den1-mtrx
          (if (or (= 0 (line-number dffr-z2-mtrx)) (= 0 (column-number dffr-z2-mtrx)))
              nil-mtrx
            dffr-z2-mtrx))         
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
     (type list fltr-p-basis  z2-list)
     (type fixnum p-basis-l den-column-n den-line-n)
     (type matrix z2-mtrx dffr-mtrx nil-mtrx 
           dffr-z2-mtrx  den1-mtrx rslt))
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

(defun mtrx-conc-2basis (mtrx1 mtrx2 basis1 basis2 cmpr)
  (declare (type matrix mtrx1 mtrx2)
           (type list basis1 basis2))
   (the matrix
     (let* ((line-n-1 (line-number mtrx1))
            (column-n-1 (column-number mtrx1))
            (column-n-2 (column-number mtrx2))
            (line-n line-n-1)
            (column-n (+ column-n-1 column-n-2)))
        (declare (type fixnum line-n-1 column-n-1 
                   column-n-2 line-n column-n))
        (let ((rslt 
                #-ACLPC
                (make-array (list line-n column-n)
                  :element-type 'fixnum
                  :initial-element 0)
                #+ACLPC
                (if (or (zerop line-n) (zerop column-n))
                   (make-array (list line-n column-n)
                     :element-type 'fixnum)
                   (make-array (list line-n column-n)
                     :element-type 'fixnum
                     :initial-element 0))))
           (declare (type matrix rslt))
           (do ((j 0 (1+ j)))
               ((>= j column-n-1))
              (declare (type fixnum j))
              (do ((i 0 (1+ i)))
                   ((>= i line-n-1))
                 (declare (type fixnum i))
                (Setf (aref rslt i j) (aref mtrx1 i j))))
          (dotimes (i (length basis2))
            (dotimes (k (length basis1))
              (if (eq :equal (funcall cmpr (nth i basis2) (nth k basis1)))

           (do ((j 0 (1+ j)))
                ((>= j column-n-2))
              (declare (type fixnum j))
              
                 (Setf (aref rslt k (+ column-n-1 j)) (aref mtrx2 i j))))))
            rslt))))


(DEFUN MULTIPRST-M-DEN-MTRX-1 (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))         
         (nil-mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (indx-list (downset-less b)))
    (declare
     (type list fltr-p-basis  indx-list)
     (type fixnum p-basis-l)
     (type matrix nil-mtrx mtrx))
    (dolist (i indx-list)
      (setf mtrx (mtrx-conc mtrx (multiprst-m-num-mtrx-1 gfltrcm p (list i) degr))))
    (let* ((den1-mtrx
            (if (or (= 0 (line-number mtrx)) (= 0 (column-number mtrx)))
                nil-mtrx
              mtrx))         
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
       (type fixnum den-column-n den-line-n)
       (type matrix den1-mtrx rslt))    
      (the matrix
        (progn
          (do ((j 0 (1+ j)))
              ((>= j den-column-n))
            (declare (type fixnum j))
            (do ((i 0 (1+ i)))
                ((>= i (min den-line-n (line-number den1-mtrx))))
              (declare (type fixnum i))
              (setf (aref rslt i j) (aref den1-mtrx i j))))
          rslt)))))


(DEFUN MULTIPRST-M-DEN-MTRX-2 (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* (
         (cmpr (cmpr gfltrcm))
         (fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         
         (nil-mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (indx-list (downset-non-comparable b t1 t2))
         (indx-list2 (append (downset-non-comparable p t1 t2) (downset-less p))))
    (declare
     (type list fltr-p-basis indx-list)
     (type fixnum p-basis-l)
     (type matrix nil-mtrx mtrx))
    (progn
      (dolist (i indx-list)
        (setf mtrx (mtrx-conc mtrx (multiprst-m-num-mtrx-1 gfltrcm p (list i) degr))))
      (dolist (i indx-list2)
        ;;(print i)
        ;;(print mtrx)
        ;;(print (multiprst-m-num-mtrx-1 gfltrcm (list i) b degr))
        (setf mtrx (mtrx-conc-2basis mtrx  (multiprst-m-num-mtrx-1 gfltrcm (list i) b degr)
                    fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr (list i)) cmpr))))
 
    (let* ((den1-mtrx
            (if (or (= 0 (line-number mtrx)) (= 0 (column-number mtrx)))
                nil-mtrx
              mtrx))         
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
       (type fixnum den-column-n den-line-n)
       (type matrix den1-mtrx rslt))    
      (the matrix
        (progn
          (do ((j 0 (1+ j)))
              ((>= j den-column-n))
            (declare (type fixnum j))
            (do ((i 0 (1+ i)))
                ((>= i (min den-line-n (line-number den1-mtrx))))
              (declare (type fixnum i))
              (setf (aref rslt i j) (aref den1-mtrx i j))))
          rslt)))))


(DEFUN MULTIPRST-M-NUM-MTRX (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))
         
         (nil-mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (mtrx (mtrx-conc (multiprst-m-num-mtrx-1 gfltrcm p b degr) 
                          (multiprst-m-den-mtrx-2 gfltrcm p b degr)))
         (den1-mtrx
          (if (or (= 0 (line-number mtrx)) (= 0 (column-number mtrx)))
              nil-mtrx
            mtrx))         
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
     (type list fltr-p-basis)
     (type fixnum p-basis-l den-column-n den-line-n)
     (type matrix nil-mtrx mtrx den1-mtrx rslt))
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


(DEFUN MULTIPRST-M-DEN-MTRX (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
         (p-basis-l (length fltr-p-basis))         
         (nil-mtrx 
          #-ACLPC (make-array (list 0 0)
                              :element-type 'fixnum
                              :initial-element 0)
          #+ACLPC (make-array (list 0 0)
                              :element-type 'fixnum))
         (mtrx (mtrx-conc (multiprst-m-den-mtrx-1 gfltrcm p b degr) 
                          (multiprst-m-den-mtrx-2 gfltrcm p b degr)))
         (den1-mtrx
          (if (or (= 0 (line-number mtrx)) (= 0 (column-number mtrx)))
              nil-mtrx
            mtrx))         
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
     (type list fltr-p-basis)
     (type fixnum p-basis-l den-column-n den-line-n)
     (type matrix nil-mtrx mtrx den1-mtrx rslt))
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


(DEFUN EFF-MULTIPRST-M-BASIS-DVS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (the list
    (let* ((fltr-p-basis (dzm-gen-fltrd-basis gfltrcm degr p))
           (num-mtrx (multiprst-m-num-mtrx gfltrcm p b degr))
           (den-mtrx (multiprst-m-den-mtrx gfltrcm p b degr)) 
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


(DEFUN EFF-MULTIPRST-M-GNRTS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((basis-dvs (multiprst-m-basis-dvs gfltrcm p b degr))
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


(DEFUN EFF-MULTIPRST-M-GROUP (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))  
  (let* ((basis-dvs (MULTIPRST-M-basis-dvs gfltrcm p b degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Multipersistence group M[~D,~D]_{~D}" p b degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))


;; USING EFFECTIVE HOMOLOGY TO DETERMINE THE NEW INVARIANT M


(DEFUN MULTIPRST-M-BASIS-DVS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-m-basis-dvs cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-multiprst-m-basis-dvs gfltrcm p b degr)
        (let ((m-list (eff-multiprst-m-basis-dvs eff-chcm p b degr)))
          (declare (list m-list))
          (if (eql gfltrcm eff-chcm) m-list
            (list
             (mapcar
                 #'(lambda (cmbn)
                     (lf hmtp-eq (rg hmtp-eq cmbn)))
               (first m-list))
             (second m-list))))))))


(DEFUN MULTIPRST-M-GNRTS (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-m-gnrts cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-multiprst-m-gnrts gfltrcm p b degr)
        (let ((gnrts (eff-multiprst-m-gnrts eff-chcm p b degr)))
          (declare (list gnrts))
          (if (eql gfltrcm eff-chcm) gnrts
            (mapcar
                #'(lambda (cmbn)
                    (lf hmtp-eq (rg hmtp-eq cmbn)))
              gnrts)))))))


(DEFUN MULTIPRST-M-GROUP (gfltrcm p b degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt p b))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "multiprst-m-group cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
        (eff-multiprst-m-group gfltrcm p b degr)
      (eff-multiprst-m-group eff-chcm p b degr))))


(DEFUN EFF-MULTIPRST-I-BASIS-DVS (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))
  (the list
    (let* ((fltr-p1-basis (gen-fltrd-basis gfltrcm degr p1))
           (num-mtrx      (multiprst-den-mtrx gfltrcm p1 p2 degr))
           (den-mtrx (mtrx-conc (multiprst-den-mtrx gfltrcm p1 s2 degr) 
                                (multiprst-den-mtrx gfltrcm s1 p2 degr))) 
           (cmbn-list nil)
           (divs nil))
      (declare
       (list fltr-p1-basis cmbn-list divs)
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
                          (mark2 fltr-p1-basis (cdr mark2)))
                         ((or (endp mark1) (endp mark2)))
                       (declare (type list mark1 mark2))
                       (if (not (= 0 (car mark1)))
                           (setq cmbn (2cmbn-add cmpr cmbn 
                                                 (cmbn degr (car mark1) (car mark2))))))
                     cmbn))) 
           cmbn-list)
         divs)))))


(DEFUN EFF-MULTIPRST-I-GNRTS (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))
  (let* ((basis-dvs (MULTIPRST-I-basis-dvs gfltrcm s1 p1 s2 p2 degr))
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


(DEFUN EFF-MULTIPRST-I-GROUP (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))  
  (let* ((basis-dvs (MULTIPRST-I-basis-dvs gfltrcm s1 p1 s2 p2 degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Multipersistence group I[~D,~D,~D,~D]_{~D}" s1 p1 s2 p2 degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))


;; USING EFFECTIVE HOMOLOGY TO DETERMINE THE NEW INVARIANT M


(DEFUN MULTIPRST-I-BASIS-DVS (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "MULTIPRST-I-basis-dvs cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-MULTIPRST-I-basis-dvs gfltrcm s1 p1 s2 p2 degr)
        (let ((m-list (eff-MULTIPRST-I-basis-dvs eff-chcm s1 p1 s2 p2 degr)))
          (declare (list m-list))
          (if (eql gfltrcm eff-chcm) m-list
            (list
             (mapcar
                 #'(lambda (cmbn)
                     (lf hmtp-eq (rg hmtp-eq cmbn)))
               (first m-list))
             (second m-list))))))))


(DEFUN MULTIPRST-I-GNRTS (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "MULTIPRST-I-gnrts cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (progn
      (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
          (eff-MULTIPRST-I-gnrts gfltrcm s1 p1 s2 p2 degr)
        (let ((gnrts (eff-MULTIPRST-I-gnrts eff-chcm s1 p1 s2 p2 degr)))
          (declare (list gnrts))
          (if (eql gfltrcm eff-chcm) gnrts
            (mapcar
                #'(lambda (cmbn)
                    (lf hmtp-eq (rg hmtp-eq cmbn)))
              gnrts)))))))


(DEFUN MULTIPRST-I-GROUP (gfltrcm s1 p1 s2 p2 degr)
  (declare 
   (type generalized-filtered-chain-complex gfltrcm)
   (type fixnum degr)
   (type gnrt s1 p1 s2 p2))
  (let* ((hmtp-eq (efhm gfltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (when (and (eq (basis gfltrcm) :locally-effective) (not (typep eff-chcm 'generalized-filtered-chain-complex)))
      (error "MULTIPRST-I-group cannot work with a LOCALLY-EFFECTIVE chain complex such that its effective chain complex is not a generalized filtered chain complex."))
    (if (not (typep eff-chcm 'generalized-filtered-chain-complex))
        (eff-MULTIPRST-I-group gfltrcm s1 p1 s2 p2 degr)
      (eff-MULTIPRST-I-group eff-chcm s1 p1 s2 p2 degr))))




