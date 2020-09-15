(in-package :cat)

(cat-init)
(setf K1
  (build-finite-ss
   '(a b c d e f
       1 ab (b a) ac (c a) bc (c b) bd (d b) cd (d c) ce (e c) de (e d) df (f d) ef (f e)
       2  bcd (cd bd bc) cde (de ce cd)
       )))

(setf k1-gflin 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a 'b 'c 'ab 'ac 'bc))
            (l12 (list 'd 'bd 'cd))
            (l13 (list 'e 'ce 'de 'bcd))
            (l21 (list 'd 'bd 'cd))
            (l22 (list 'e 'ce 'de))
            (l23 (list 'cde))
            (l31 (list 'e 'ce 'de))
            (l32 (list 'bcd 'cde))
            (l33 (list 'f 'df 'ef))
            (rslt +empty-list+))
        (declare (list l11 l12 l13 l21 l22 l23 l31 l32 l33 rslt))
        (progn
          (if (contains l11 gnrt)  (push '(1 1) rslt))
          (if (contains l12 gnrt)  (push '(1 2) rslt))
          (if (contains l13 gnrt)  (push '(1 3) rslt))
          (if (contains l21 gnrt)  (push  '(2 1) rslt))     
          (if (contains l22 gnrt) (push  '(2 2) rslt))
          (if (contains l23 gnrt) (push  '(2 3) rslt))
          (if (contains l31 gnrt) (push '(3 1) rslt))
          (if (contains l32 gnrt) (push  '(3 2) rslt))
          (if (contains l33 gnrt) (push  '(3 3) rslt))
          (nreverse rslt)))))


(setf K (change-chcm-to-gflcc K1 (z2) k1-gflin 'k1-gflin))

(multiprst-group K '(1 2) '(2 2) 1)

(multiprst-gnrts K '(1 2) '(2 2) 1)

(multiprst-group K '(1 2) '(3 3) 1)

(multiprst-gnrts K '(1 2) '(3 3) 1)



;; We build the second example and we test that the rank invariant of filtrations in Figures 1 & 2 are the same


(setf K21
  (build-finite-ss
   '(a b c d e f
       1 ab (b a) ac (c a) bc (c b) bd (d b) cd (d c) ce (e c) de (e d) df (f d) ef (f e)
       2  bcd (cd bd bc) cde (de ce cd)
       )))

(setf k21-gflin 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a 'b 'c 'ab 'ac 'bc))
            (l12 (list 'd 'bd 'cd))
            (l13 (list 'e 'ce 'de 'bcd))
            (l21 (list 'd 'e 'cd 'de 'ce))
            (l22 (list 'bd))
            (l23 (list 'cde))
            (l31 (list 'bd))
            (l32 (list 'bcd 'cde))
            (l33 (list 'f 'df 'ef))
            (rslt +empty-list+))
        (declare (list l11 l12 l13 l21 l22 l23 l31 l32 l33 rslt))
        (progn
          (if (contains l11 gnrt)  (push '(1 1) rslt))
          (if (contains l12 gnrt)  (push '(1 2) rslt))
          (if (contains l13 gnrt)  (push '(1 3) rslt))
          (if (contains l21 gnrt)  (push  '(2 1) rslt))     
          (if (contains l22 gnrt) (push  '(2 2) rslt))
          (if (contains l23 gnrt) (push  '(2 3) rslt))
          (if (contains l31 gnrt) (push '(3 1) rslt))
          (if (contains l32 gnrt) (push  '(3 2) rslt))
          (if (contains l33 gnrt) (push  '(3 3) rslt))
          (nreverse rslt)))))

(setf K2 (change-chcm-to-gflcc K21 (z2) k21-gflin 'k21-gflin))

(multiprst-group K '(1 2) '(2 2) 1)

(multiprst-gnrts K2 '(1 2) '(2 2) 1)

(multiprst-group K2 '(1 2) '(3 3) 1)

(multiprst-gnrts K2 '(1 2) '(3 3) 1)

(multiprst-gnrts K '(1 2) '(3 2) 1))

(multiprst-gnrts K2 '(1 2) '(3 2) 1)

(length (multiprst-gnrts K '(1 2) '(3 2) 1))

(length (multiprst-gnrts K2 '(1 2) '(3 2) 1))

(multiprst-group K (list 1 1) (list 1 1) 1)


;; We compute the new descriptor and the new invariant


(setf k1-gflin2 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a 'b 'c 'ab 'ac 'bc))
            (l12 (list 'd 'bd 'cd))
            (l13 (list 'e 'ce 'de 'bcd))
            (l21 (list 'd 'bd 'cd))
            (l22 (list 'e 'ce 'de))
            (l23 (list 'cde))
            (l31 (list 'e 'ce 'de))
            (l32 (list 'bcd 'cde))
            (l33 (list 'f 'df 'ef))
            (rslt +empty-list+))
        (declare (list l11 l12 l13 l21 l22 l23 l31 l32 l33 rslt))
        (progn
          (if (contains l11 gnrt)  (push (list '(1 1)) rslt))
          (if (contains l12 gnrt)  (push (list '(1 2)) rslt))
          (if (contains l13 gnrt)  (push (list '(1 3)) rslt))
          (if (contains l21 gnrt)  (push  (list '(2 1)) rslt))     
          (if (contains l22 gnrt) (push  (list '(2 2)) rslt))
          (if (contains l23 gnrt) (push  (list '(2 3)) rslt))
          (if (contains l31 gnrt) (push (list '(3 1)) rslt))
          (if (contains l32 gnrt) (push (list  '(3 2)) rslt))
          (if (contains l33 gnrt) (push (list  '(3 3)) rslt))
          (nreverse rslt)))))

(setf K3 (change-chcm-to-gflcc K1 (downsets (z2)) k1-gflin2 'k1-gflin2))
(setf t1 3 t2 3)


(setf k21-gflin2 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a 'b 'c 'ab 'ac 'bc))
            (l12 (list 'd 'bd 'cd))
            (l13 (list 'e 'ce 'de 'bcd))
            (l21 (list 'd 'e 'cd 'de 'ce))
            (l22 (list 'bd))
            (l23 (list 'cde))
            (l31 (list 'bd))
            (l32 (list 'bcd 'cde))
            (l33 (list 'f 'df 'ef))
            (rslt +empty-list+))
        (declare (list l11 l12 l13 l21 l22 l23 l31 l32 l33 rslt))
        (progn
          (if (contains l11 gnrt)  (push (list '(1 1)) rslt))
          (if (contains l12 gnrt)  (push (list '(1 2)) rslt))
          (if (contains l13 gnrt)  (push (list '(1 3)) rslt))
          (if (contains l21 gnrt)  (push  (list '(2 1)) rslt))     
          (if (contains l22 gnrt) (push  (list '(2 2)) rslt))
          (if (contains l23 gnrt) (push  (list '(2 3)) rslt))
          (if (contains l31 gnrt) (push (list '(3 1)) rslt))
          (if (contains l32 gnrt) (push  (list '(3 2)) rslt))
          (if (contains l33 gnrt) (push  (list '(3 3)) rslt))
          (nreverse rslt)))))

(setf K4 (change-chcm-to-gflcc K1 (downsets (z2)) k21-gflin2 'k21-gflin2))

(multiprst-m-group K3 (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)
(multiprst-m-gnrts K3 (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)

(multiprst-m-group K4 (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)

(multiprst-m-group K3 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)

(multiprst-m-gnrts K3 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)
(multiprst-m-group K4 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)

(multiprst-m-group K4 (list '(1 3) '(2 1) ) (list '(2 3) '(3 2)) 1)
(multiprst-m-gnrts K4 (list '(1 3) '(2 1) ) (list '(2 3) '(3 2)) 1)

(multiprst-m-group K4 (list '(1 3) '(2 1)) (list '(2 3) '(3 2)) 1)
(multiprst-m-group K3 (list '(1 3) '(2 1)) (list '(2 3) '(3 2)) 1)

(multiprst-i-group K3 (list '(1 1) ) (list '(1 2) '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)
(multiprst-i-group K4 (list '(1 1) ) (list '(1 2) '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)

(multiprst-i-gnrts K3 (list '(1 1) ) (list '(1 2) '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)
(multiprst-i-gnrts K4 (list '(1 1) ) (list '(1 2) '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)

(multiprst-i-group K3 (list '(1 1) ) (list  '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)
(multiprst-i-group K4 (list '(1 1) ) (list '(2 1)) (list '(1 2) '(2 1)) (list '(1 3) '(3 2)) 1)

(multiprst-i-group K3 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)

(multiprst-i-gnrts K3 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)
(multiprst-i-group K4 (list '(1 3) '(2 2) '(3 1)) (list '(2 3) '(3 2)) 1)

(multiprst-i-group K4 (list '(1 3) '(2 1) ) (list '(2 3) '(3 2)) 1)
(multiprst-i-gnrts K4 (list '(1 3) '(2 1) ) (list '(2 3) '(3 2)) 1)

(multiprst-i-group K4 (list '(1 3) '(2 1)) (list '(2 3) '(3 2)) 1)
(multiprst-i-group K3 (list '(1 3) '(2 1)) (list '(2 3) '(3 2)) 1)

(multiprst-i-group K3  (list '(1 2) '(2 1)) (list '(2 2)) (list '(2 3) '(3 2)) (list '(3 3)) 1)
(multiprst-i-group K4  (list '(1 2) '(2 1)) (list '(2 2)) (list '(2 3) '(3 2)) (list '(3 3)) 1)



;; Example of filtered chain complex (not corresponding to a simplicial set)

(cat-init)
(setf chcm
  (the chain-complex
    (build-chcm
     :cmpr #'s-cmpr
     :basis #'(lambda (dmns)
                (the list
                  (case dmns (0 (list 'a)) 
                    (1 (list 'b1 'b2 'b3 'b4 ))
                    (2 (list 'c1 'c2 'c3 'c4))
                    (otherwise +empty-list+))))
     :bsgn 'a
     :intr-dffr #'(lambda (dmns gnrt)
                    (let ((rslt 
                           (if (and (= dmns 1) (eq gnrt 'b1)) (cmbn 0 1 'a)
                             (if (and (= dmns 2)  (eq gnrt 'c1)) (cmbn 1 2 'b2 1 'b4)
                               (if (and (= dmns 2)  (eq gnrt 'c2)) (cmbn 1 2 'b3)
                                 (if (and (= dmns 2)  (eq gnrt 'c3)) (cmbn 1 2 'b4)
                                   (if (and (= dmns 2)  (eq gnrt 'c4)) (cmbn 1 1 'b3))))))))
                      (if rslt rslt (cmbn (1- dmns)))))
                                 
     :strt :gnrt
     :orgn '(chcm1))))


(setf chcm-gflin 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a))
            (l12 (list 'b1 'b2 'b3 'b4 'c1 'c2 'c3))
            (l21 (list 'b3 'b4  'c3))
            (l22 (list 'c4))
            (rslt +empty-list+))
        (declare (list l11 l12 l21 l11 rslt))
        (progn
          (if (contains l11 gnrt)  (push '(1 1) rslt))
          (if (contains l12 gnrt)  (push '(1 2) rslt))
          (if (contains l21 gnrt)  (push '(2 1) rslt))
          (if (contains l22 gnrt)  (push '(2 2) rslt))
          (nreverse rslt)))))

(setf C (change-chcm-to-gflcc chcm (z2) chcm-gflin 'chcm-gflin))

(multiprst-group C '(1 1) '(1 1) 0)
;;Multipersistence group H[(1 1),(1 1)]_{0}
;;Component Z
(multiprst-group C '(1 1) '(2 2) 0)
;;Multipersistence group H[(1 1),(2 2)]_{0}
;;NIL
(multiprst-group C '(1 1) '(1 2) 0)
;;Multipersistence group H[(1 1),(1 2)]_{0}
;;NIL
(multiprst-group C '(1 1) '(2 1) 0)
;;Multipersistence group H[(1 1),(2 1)]_{0}
;;Component Z
(multiprst-group C '(1 1) '(1 1) 1)
;;Multipersistence group H[(1 1),(1 1)]_{1}
;;NIL
(multiprst-group C '(1 2) '(1 2) 1)
;;Multipersistence group H[(1 2),(1 2)]_{1}
;;Component Z/2Z
;;Component Z/4Z
(multiprst-group C '(2 1) '(2 1) 1)
;;Multipersistence group H[(2 1),(2 1)]_{1}
;;Component Z/2Z
;;Component Z
(multiprst-group C '(2 1) '(2 2) 1)
;;Multipersistence group H[(2 1),(1 1)]_{1}
;;Component Z/2Z
(multiprst-group C '(2 2) '(2 2) 1)
;;Multipersistence group H[(2 2),(2 2)]_{1}
;;Component Z/4Z
(multiprst-group C '(2 2) '(2 2) 2)
;;Multipersistence group H[(1 1),(1 1)]_{2}
;;Component Z


(multiprst-gnrts C '(1 2) '(1 2) 1)
(multiprst-gnrts C '(1 0) '(1 0) 1)
(multiprst-gnrts C '(1 1) '(1 1) 1)
(multiprst-gnrts C '(1 0) '(1 1) 1)


;; We compute the new descriptor and invariant

(setf chcm-gflin2 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l11 (list 'a))
            (l12 (list 'b1 'b2 'b3 'b4 'c1 'c2 'c3))
            (l21 (list 'b3 'b4  'c3))
            (l22 (list 'c4))
            (rslt +empty-list+))
        (declare (list l11 l12 l21 l11 rslt))
        (progn
          (if (contains l11 gnrt)  (push (list '(1 1)) rslt))
          (if (contains l12 gnrt)  (push (list '(1 2)) rslt))
          (if (contains l21 gnrt)  (push (list '(2 1)) rslt))
          (if (contains l22 gnrt)  (push (list '(2 2)) rslt))
          (nreverse rslt)))))


(setf C2 (change-chcm-to-gflcc chcm (downsets (z2)) chcm-gflin2 'chcm-gflin3))
(setf C C2)

(setf t1 2 t2 2)
(multiprst-m-group C (list '(1 1)) (list '(2 1)) 0)
(multiprst-m-group C (list '(1 1)) (list '(1 2)) 0)

(multiprst-m-group C (list '(1 2)) (list '(2 2)) 1)
(multiprst-m-group C (list '(2 1)) (list '(2 2)) 1)
(multiprst-m-group C (list '(1 2) '(2 1)) (list '(2 2)) 1)
(multiprst-m-gnrts C (list '(1 2) '(2 1)) (list '(2 2)) 1)

(multiprst-i-group C (list '(1 1)) (list '(1 2)) (list '(1 2)) (list '(2 2)) 1)
(multiprst-i-group C (list '(1 1)) (list '(2 1)) (list '(2 1)) (list '(2 2)) 1)
(multiprst-i-group C (list '(1 1)) (list '(1 2) '(2 1)) (list '(1 2) '(2 1)) (list '(2 2)) 1)


;; Tower of Postnikov of the 3-sphere
(cat-init)
(setf s3 (sphere 3))
(setf k3 (chml-clss s3 3))
(setf F3 (z-whitehead s3 k3))
(setf X4 (fibration-total F3))
(setf k4 (chml-clss x4 4))
(setf F4 (z2-whitehead X4 k4))
(setf X5 (fibration-total F4))
(setf k5 (chml-clss x5 5))
(setf F5 (z2-whitehead X5 k5))
(setf X6aux (fibration-total F5))
(setf ex6 (rbcc (efhm X6aux)))
(setf X6 (change-chcm-to-gflcc ex6 (dzn 3) 3tnpr-gflin '3tnpr-gflin))


(multiprst-group X6 (list '(7 7 7)) (list '(7 7 7)) 6)

(multiprst-group X6 (list '(3 0 0)) (list '(7 7 7)) 6)

(multiprst-group X6 (list '(0 0 3)) (list '(0 0 3)) 6)

(multiprst-group X6 (list '(3 0 0)) (list '(3 0 0)) 6)

(multiprst-group X6 (list '(0 0 3)) (list '(0 0 3)) 3)



;; Using discrete vector fields on digital images

(cat-init)
(setf hlist nil vlis nil)

(dotimes (i 4)
      (format t "Select the horizontal image ~D~%:" 
        i)
  (push (image-to-ss) hlist))
(setf hlist (nreverse hlist))

(setf vlist nil)
(dotimes (i 4)
      (format t "Select the vertical image ~D~%:" 
        i)
  (push (image-to-ss) vlist))
(setf vlist (nreverse vlist))


  (setf l nil)
(setf l (list-of-list hlist vlist))

(cat-init)

(setf ss (gfltrcm-from-list l))

(setf rslts-list (build-empty-gen-rslts))
(setf (slot-value ss 'gen-rslts) rslts-list)

(setf gfltrcm ss mtrx-degr 1)
(setf vf-list (gfltrcm-vf-list gfltrcm mtrx-degr))

(setf crtc-basis-0 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr vf-list 0))
(setf crtc-basis-1 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr vf-list 1))
(setf crtc-basis-2 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr vf-list 2))

(setf crtc-basis #'(lambda (degr)
                     (if (eq 0 degr) crtc-basis-0
                       (if (eq 1 degr) crtc-basis-1
                         (if (eq 2 degr) crtc-basis-2
                           nil)))))

(setf r42 (gfltr-vf-1degr-rdct2 ss 1))
(setf crtc42 (bcc r42))

(setf gfltrcm crtc42 mtrx-degr2 2)
(setf vf-list2 (gfltrcm-vf-list gfltrcm mtrx-degr2))

(setf crtc-basis2-0 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr2 vf-list2 0))
(setf crtc-basis2-1 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr2 vf-list2 1))
(setf crtc-basis2-2 (gfltrcm-vf-crtc-list gfltrcm mtrx-degr2 vf-list2 2))

(setf crtc-basis2 #'(lambda (degr)
                     (if (eq 0 degr) crtc-basis2-0
                       (if (eq 1 degr) crtc-basis2-1
                         (if (eq 2 degr) crtc-basis2-2
                           nil)))))

(setf r43 (gfltr-vf-1degr-rdct3 gfltrcm 2))
(setf crtc43 (bcc r43))

(setf efchcm (build-gflcc 
              :cmpr (cmpr crtc43)
              :basis (basis crtc43)
              :bsgn (bsgn crtc43)
              :intr-dffr #'(lambda (degr gnrt)
                             (? (dffr crtc43) degr gnrt))
              :dffr-strt :gnrt 
              :poset (z2)
              :gen-flin (gen-flin crtc43)
              :orgn 'efchcm))

(length (basis efchcm 0))
(length (basis efchcm 1))
(length (basis efchcm 2))

(homology efchcm 1)

(multiprst-group efchcm '(1 1) '(2 2) 1)
(multiprst-group efchcm '(1 1) '(3 3) 1)

(setf K5 (build-gflcc 
              :cmpr (cmpr crtc43)
              :basis (basis crtc43)
              :bsgn (bsgn crtc43)
              :intr-dffr #'(lambda (degr gnrt)
                             (? (dffr crtc43) degr gnrt))
              :dffr-strt :gnrt 
              :poset (downsets (z2))
               :gen-flin #'(lambda (degr gnrt)
                             (list (funcall (gen-flin crtc43) degr gnrt)))
          :orgn 'efchcm2))

(multiprst-m-group K5 (list '(1 0)) (list '(1 3) '(3 0)) 1)
(multiprst-m-group K5 (list '(0 2) '(1 0)) (list '(0 3) '(3 0)) 1)

(multiprst-group K5 (list  '(0 0)) (list '(0 3) '(3 0)) 1)
(multiprst-group K5 (list  '(0 1)) (list '(0 1)) 1)




