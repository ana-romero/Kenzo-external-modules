;; GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF
;; GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF
;; GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF   GENERALIZED-FILTERED-VF


(DEFUN IS-FIRST-IN-PAIR (n pair)
  (declare (fixnum n)
           (list pair))
  (eq n (first pair)))


(DEFUN IS-SECOND-IN-PAIR (n pair)
  (declare (fixnum n)
           (list pair))
  (eq n (second pair)))


(DEFUN MATRICE-TO-LMTRX (mtrx)
  (let* (
	  (column-n (ncol mtrx))
        (rslt (make-array column-n)))
    (dotimes (j column-n)
      (let ((rsltj '())
            (columnj-pair-list 
             (let ((ptc (basecol mtrx (1+ j)))
                   (res '()))
               (do ((pc (up ptc) (up pc)))
                   ((eq pc ptc))
                 (push (list (ilig pc) (val pc)) res))
               res)))
        
        (mapcar #'(lambda (pair)
                    (push  (list (1- (first pair)) (second pair)) rsltj))
         columnj-pair-list)
         
        (setf (svref rslt j) (nreverse rsltj) )))
    rslt))


(DEFUN GFLTR-INDEX-LIST (gfltrcm degr)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (fixnum degr))
  (the list
    (let ((rslt nil))
      (declare 
       (type list rslt))
      (mapcar #'(lambda (gnrt)
                  (declare (type gnrt gnrt))
                  (let ((gflin (gen-flin gfltrcm degr gnrt)))
                    (declare (type list gflin))
                    (if (not (position gflin rslt :test #'(lambda (l1 l2)
                                                            (eq :equal (pocmpr (downsets (pos gfltrcm)) l1 l2)))))
                        (push gflin rslt))))       
        (basis gfltrcm degr))
      rslt)))


(DEFUN GFLCC-DFFR-SUBMATRICE (gfltrcm degr gen-flin)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (fixnum degr)
   (list gen-flin))
  (when (eq (basis gfltrcm) :locally-effective)
    (error "GFLCC-DFFR-SUBMATRICE cannot work with a LOCALLY-EFFECTIVE chain complex."))
  (the matrice
    (let* ((cmpr (cmpr1 gfltrcm))
           (sbasis (gen-fltrd-gflin-elements gfltrcm degr gen-flin))
           (tbasis (gen-fltrd-gflin-elements gfltrcm  (1- degr) gen-flin))
           (srank (length sbasis))
           (trank (length tbasis)))
      (declare
       (type cmprf cmpr)
       (list sbasis tbasis)
       (fixnum srank trank))
      (let (
            (mat (creer-matrice trank srank))
            (test #'(lambda (gnrt1 gnrt2)
                      (eq (funcall cmpr gnrt1 gnrt2) :equal))))
        (declare
         
         (type matrice mat)
         (function test))
        (do ((i 1 (1+ i))
             (mark sbasis (cdr mark)))
            ((endp mark))
          (declare
           (fixnum i)
           (list mark))
          (maj-colonne mat i
                       (mapcar #'(lambda (term)
                                   (declare (type term term))
                                   (list
                                    (1+ (position (gnrt term) tbasis :test test))
                                    (cffc term)))
                         (let ((l (cmbn-list (? gfltrcm degr (car mark))))
                               (rslt nil))
                           (mapcar #'(lambda (term)
                                       (declare (type term term))
                                       (if (find (gnrt term) tbasis :test test)
                                           (push term rslt)))
                             l)
                           rslt))))
        mat))))


(DEFUN GFLTRCM-VF-LIST (gfltrcm degr)
  (let ((indx-list (gfltr-index-list gfltrcm degr)))
    (mapcar #'(lambda (gflin)
                (list gflin
                      (m-vf-hard (length (gen-fltrd-gflin-elements gfltrcm  (1- degr) gflin)) 
                                 (matrice-to-lmtrx (GFLCC-DFFR-SUBMATRICE gfltrcm degr gflin)))))
      indx-list)))


(DEFUN GFLTRCM-VF (gfltrcm mtrx-degr vf-list)
  (flet ((rslt (dmns gmsm)
               (if (or (> dmns mtrx-degr) (< dmns (1- mtrx-degr)))
                   (vctr :crtc)
                 (let* ((gflin (gen-flin gfltrcm mtrx-degr gmsm))
                        (lmtrx (matrice-to-lmtrx (GFLCC-DFFR-SUBMATRICE gfltrcm mtrx-degr gflin)))
                        (gflin-vf-list (second (find gflin vf-list :test #'(lambda (l1 l2)
                                                                             (eq :equal (pocmpr (downsets (pos gfltrcm)) l1 (first l2)))))))
                        (gflin-elts-degr (GEN-FLTRD-GFLIN-ELEMENTS gfltrcm mtrx-degr gflin))
                        (gflin-elts-degr-1 (GEN-FLTRD-GFLIN-ELEMENTS gfltrcm (1- mtrx-degr) gflin))
                        (cmpr 
                         #'(lambda (gnrt1 gnrt2)
                             (eql :equal (cmpr gfltrcm gnrt1 gnrt2))))
                        )
                   (if (eq (1- mtrx-degr) dmns)
                       (let* ((i (position gmsm gflin-elts-degr-1 :test cmpr))
                              (p (find i gflin-vf-list :test #'IS-FIRST-IN-PAIR)))
                         (if p 
                             (let* ((im-p (second p))
                                    (lmtrx-column-im-p (svref  lmtrx im-p ))
                                    (img (second (find i lmtrx-column-im-p :test #'IS-FIRST-IN-PAIR))))
                               (vctr :sorc  (nth im-p gflin-elts-degr) img))
                           (vctr :crtc)))
                     (let* ((i (position gmsm gflin-elts-degr :test cmpr))
                            (p (find i gflin-vf-list :test #'IS-SECOND-IN-PAIR)))
                       (if p 
                           (let* ((im-p (first p))
                                  (lmtrx-column-i (svref lmtrx i))
                                  (img (second (find im-p lmtrx-column-i :test #'IS-FIRST-IN-PAIR))))
                             (vctr :trgt (nth im-p gflin-elts-degr-1) img))
                         (vctr :crtc))))))))
    #'rslt))


(DEFUN GFLTRCM-VF-CRTC-LIST (gfltrcm mtrx-degr vf-list degr)
  (if (not (or (eq degr mtrx-degr) (eq degr (1- mtrx-degr))))
      (basis gfltrcm degr)
    (let ((rslt +empty-list+)
          (indx-list (GFLTR-INDEX-LIST gfltrcm degr))
          (cmpr (if (eq (1- mtrx-degr) degr) #'IS-FIRST-IN-PAIR
                  #'IS-SECOND-IN-PAIR)))
      (progn
        (mapcar #'(lambda (gflin)
                    (let ((gflin-elts-degr (GEN-FLTRD-GFLIN-ELEMENTS gfltrcm degr gflin))
                          (gflin-vf-list (second (find gflin vf-list :test #'(lambda (l1 l2)
                                                                               (eq :equal (pocmpr (downsets (pos gfltrcm)) l1 (first l2))))))))
                      (dotimes (i (length gflin-elts-degr))
                        (if (not (find i gflin-vf-list :test cmpr)) (push (nth i gflin-elts-degr) rslt)))))
          indx-list)
        
        (nreverse rslt)))))



(DEFUN GFLTRCM-VF-IS-CRTC (gfltrcm mtrx-degr vf-list degr gnrt)
  (let ((crtc-list (gfltrcm-vf-crtc-list gfltrcm mtrx-degr vf-list degr))
        (cmpr 
         #'(lambda (gnrt1 gnrt2)
             (eql :equal (cmpr gfltrcm gnrt1 gnrt2)))))
    (if (find gnrt crtc-list :test cmpr) 't nil)))


(DEFUN GFLTRCM-CRTC-BASIS (gfltrcm mtrx-degr vf-list)
  (flet ((rslt (degr)
               (gfltrcm-vf-crtc-list gfltrcm mtrx-degr vf-list degr)))
    (the basis #'rslt)))


(DEFUN GFLTRCM-CRTC-INTR-DFFR (chcm vf &aux
                                      (dffr (dffr chcm))
                                      (hh (chcm-vf-reduction-h chcm vf))
                                      (cmpr (cmpr chcm)))
  (declare (type chain-complex chcm)
           (type vector-field vf)
           (type morphism dffr hh)
           (type cmprf cmpr))
  (the intr-mrph
    (flet
        ((rslt
          (dmns gmsm)
                    
          (the cmbn
               (2cmbn-sbtr
               cmpr
               (cmbn-select-stts vf
                                    (? dffr dmns gmsm)
                                    :crtc)
                (cmbn-select-stts vf
                      (? dffr (? hh
                           (cmbn-select-stts vf
                                    (? dffr dmns gmsm)
                                    :sorc)) ) :crtc)))))
      (declare (ftype (function (cmbn) cmbn) rslt))
      #'rslt)))


(DEFUN FILTERED-VF-CRTC-CHCM (gfltrcm mtrx-degr vf-list)
  (build-gflcc
   :cmpr (cmpr gfltrcm)
   :basis (gfltrcm-crtc-basis gfltrcm mtrx-degr vf-list)
   :bsgn (bspn gfltrcm)
   :intr-dffr (gfltrcm-crtc-intr-dffr gfltrcm (gfltrcm-vf gfltrcm mtrx-degr vf-list))
   :dffr-strt :gnrt
   :poset (pos gfltrcm)
   :gen-flin (gen-flin gfltrcm)
   :orgn  `(critical-chcm5 of ,gfltrcm)))


(DEFUN FILTERED-VF-ISOG (cmbn) cmbn)
(DEFUN FILTERED-VF-ISOF (cmbn) cmbn)


(DEFUN GFLTR-VF-1DEGR-RDCT (gfltrcm degr)
  (let* ((vf-list (gfltrcm-vf-list gfltrcm degr))
         (crtc-chcm (filtered-vf-crtc-chcm gfltrcm degr vf-list))
         (vf (gfltrcm-vf gfltrcm degr vf-list)))
    (chcm-vf-reduction gfltrcm vf crtc-chcm #'filtered-vf-ISOF #'filtered-vf-ISOG)))


(DEFUN GFLTR-VF-RDCT (gfltrcm degr)
  (let* ((r (gfltr-vf-1degr-rdct gfltrcm 1))
         (crtc-chcm (bcc r)))
    (dotimes (i (1- degr))
      (progn
        (setq r (cmps  (gfltr-vf-1degr-rdct crtc-chcm (+ 2 i)) r))
        (setq crtc-chcm (bcc r))))
    r))



