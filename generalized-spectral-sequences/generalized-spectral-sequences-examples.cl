(in-package :user)
(load "cat-init")
(load-cfiles)
(load "filtered-complexes.cl")
(load "spectral-sequences.cl")
(load "generalized-spectral-sequences.cl")



;;; DEFINITION OF A SIMPLICIAL SET

(setf ss
  (build-finite-ss
   '(a b c d e
     1 ab (b a) ac (c a) bc (c b) bd (d b) cd (d c) be (e b) de (e d)
        )))
;; [K1 Simplicial-Set]

;;; DEFINITION OF A GENERALIZED FILTRATION OVER Z2

(setf ss-gflin1 
  #'(lambda (degr gnrt)
      (declare 
       (type fixnum degr)
       (type gnrt gnrt))
      (let ((l00 (list 'a 'b 'c 'ab 'ac 'bc))
            (l01 (list 'd 'e 'bd 'be 'de))
            (l10 (list 'd 'bd 'cd))
            (rslt +empty-list+))
        (declare (list l00 l01 l10 rslt))
        (progn
          (if (contains l00 gnrt)  (push '(0 0) rslt))
          (if (contains l01 gnrt)  (push '(0 1) rslt))
          (if (contains l10 gnrt)  (push '(1 0) rslt))
          (nreverse rslt)))))
;; #<Interpreted Function (unnamed) @ #x20f8fed2>            

(setf ss1 (change-chcm-to-gflcc ss (z2) ss-gflin1 'ss-gflin1))
;; [K7 Generalized-Filtered-Chain-Complex]


;;; EXAMPLES OF COMPUTATIONS OF THE GENERALIZED FILTRATION OVER THE SIMPLICES OF ss1


(gen-flin ss1 0 'a)
;; ((0 0))

(gen-flin ss1 0 'd)
;; ((0 1) (1 0))

(gen-flin ss1 1 'bd)
;; ((0 1) (1 0))

(gen-flin ss1 1 'cd)
;; ((1 0))


;;; EXAMPLES OF THE GENERALIZED FILTERED BASIS OF ss1 FOR DIFFERENT FILTRATION INDEXES

(gen-fltrd-basis ss1 0 '(0 0))
;; (A B C)

(gen-fltrd-basis ss1 0 '(0 1))
;; (A B C D E)

(gen-fltrd-basis ss1 1 '(1 0))
;; (AB AC BC BD CD)

(gen-fltrd-basis ss1 1 '(1 1)) 
;; (AB AC BC BD BE CD DE)


;;; TWO EXAMPLES OF DEFINITION OF GENERALIZED FILTRATION OVER DZ2

(setf ss-gflin2 
      #'(lambda (degr gnrt)
          (declare 
           (type fixnum degr)
           (type gnrt gnrt))
          (let ((l00 (list 'a 'b 'c 'ab 'ac 'bc))
                (l01 (list 'd 'e 'bd 'be 'de))
                (l10 (list 'd 'bd 'cd))
                (rslt +empty-list+))
            (declare (list l00 l01 l10 rslt))
            (progn
              (if (contains l00 gnrt)  (push (list '(0 0)) rslt))
              (if (contains l01 gnrt)  (push (list '(0 1)) rslt))
              (if (contains l10 gnrt)  (push (list '(1 0)) rslt))
              (nreverse rslt)))))
;; #<Interpreted Function (unnamed) @ #x20e8bdba>

(setf ss2 (change-chcm-to-gflcc ss (dz2) ss-gflin2 'ss-gflin2))
;; [K10 Generalized-Filtered-Chain-Complex]


(setf ss-gflin3 
      #'(lambda (degr gnrt)
          (declare 
           (type fixnum degr)
           (type gnrt gnrt))
          (let ((l00 (list 'a 'b 'c 'ab 'ac 'bc))
                (l01 (list 'd 'bd 'cd))
                (l10 (list 'd 'bd 'cd))
                (l11 (list 'e 'be 'de))
                (rslt +empty-list+))
            (declare (list l00 l01 l10 l11 rslt))
            (progn
              (if (contains l00 gnrt)  (push (list '(0 0)) rslt))
              (if (contains l01 gnrt)  (push (list '(0 1)) rslt))
              (if (contains l10 gnrt)  (push (list '(1 0)) rslt))
              (if (contains l11 gnrt)  (push (list '(1 1)) rslt))
              (nreverse rslt)))))
;; #<Interpreted Function (unnamed) @ #x21093dea>

(setf ss3 (change-chcm-to-gflcc ss (dz2) ss-gflin3 'ss-gflin3))
;;[K12 Generalized-Filtered-Chain-Complex]



;;; EXAMPLES OF THE GENERALIZED FILTRATION OVER THE SIMPLICES OF ss2 AND ss3
;;; IN THOSE CASE THE FILTRATION INDEXES ARE DOWNSETS GIVEN BY A LIST OF ELEMENTS OF Z2


(gen-flin ss2 0 'a)
;; (((0 0)))

(gen-flin ss3 0 'd)
;; (((0 1)) ((1 0)))

(gen-flin ss2 1 'bd)
;; (((0 1)) ((1 0)))

(gen-flin ss3 1 'cd)
;; (((0 1)) ((1 0)))


;;; EXAMPLES OF THE GENERALIZED FILTERED BASIS OF ss2 AND ss3 FOR DIFFERENT FILTRATION INDEXES

(gen-fltrd-basis ss2 0 '((0 0)))
;; (A B C)

(gen-fltrd-basis ss3 0 '((0 1) (1 0)))
;; (A B C D)

(gen-fltrd-basis ss2 1 '((1 0)))
;; (AB AC BC BD CD)

(gen-fltrd-basis ss2 1 '((0 1) (1 0)))
;; (AB AC BC BD BE CD DE)

(gen-fltrd-basis ss3 1 '((1 1))) 
;; (AB AC BC BD BE CD DE)


;; COMPUTATION OF THE GENERALIZED SPECTRAL SEQUENCES OF ss1, ss2 AND ss3

(gen-spsq-group ss1 '(0 0) '(0 0) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(0 0),(0 0),(1 1),(1 1)]_{1}
;; Component Z
;; Component Z

(gen-spsq-group ss1 '(0 0) '(0 1) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(0 0),(0 1),(1 1),(1 1)]_{1}
;; Component Z

(gen-spsq-group ss1 '(1 0) '(1 0) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(1 0),(1 0),(1 1),(1 1)]_{1}
;; Component Z

(gen-spsq-group ss1 '(1 -1) '(1 -1) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(1 -1),(1 -1),(1 1),(1 1)]_{1}
;; Component Z
;; Component Z
;; Component Z

(gen-spsq-group ss1 '(0 0) '(0 0) '(1 1) '(1 1) 2)
;; Generalized spectral sequence S[(0 0),(0 0),(1 1),(1 1)]_{2}
;; NIL


(gen-spsq-group ss2 (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((1 0) (0 1)),((1 1)),((1 1))]_{1}
;; NIL

(gen-spsq-group ss2 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 1)),((1 1))]_{1}
;; Component Z
;; Component Z

(gen-spsq-group ss2 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 0) (0 1)),((1 1))]_{1}
;; Component Z
;; Component Z

(gen-spsq-group ss3 (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((1 0) (0 1)),((1 1)),((1 1))]_{1}
;; Component Z

(gen-spsq-group ss3 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
;;Generalized spectral sequence S[((0 0)),((0 0)),((1 1)),((1 1))]_{1}
;; Component Z
;; Component Z

(gen-spsq-group ss3 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 0) (0 1)),((1 1))]_{1}
;; Component Z


(gen-spsq-gnrts ss1 '(1 -1) '(1 -1) '(1 1) '(1 1) 1)
;; (
;; ----------------------------------------------------------------------{CMBN 1}
;; <1 * AB>
;; <-1 * AC>
;; <1 * BC>
;; ------------------------------------------------------------------------------
;;
;; ----------------------------------------------------------------------{CMBN 1}
;; <-1 * AB>
;; <1 * AC>
;; <-1 * BD>
;; <1 * CD>
;; ------------------------------------------------------------------------------
;; 
;; ----------------------------------------------------------------------{CMBN 1}
;; <1 * BD>
;; <-1 * BE>
;; <1 * DE>
;; ------------------------------------------------------------------------------
;; )

(gen-spsq-dffr ss1 '(0 0) '(0 1) '(1 1) '(1 1) '(-1 0) '(-1 0) '(1 1) '(1 1) 1 '(1))
;; (0)


;;; ANOTHER EXAMPLE OF SIMPLICIAL SET WITH A GENERALIZED FILTRATION OVER DZ2

(setf ss4
  (build-finite-ss
   '(a b c d e
       1 ab (b a) ac (c a) bc (c b) bd (d b) cd (d c) be (e b) de (e d)
       2 abc (bc ac ab) bcd (cd bd bc) bde (de be bd)
        )))
;; [K26 Simplicial-Set]


(setf ss-gflin4 
      #'(lambda (degr gnrt)
          (declare 
           (type fixnum degr)
           (type gnrt gnrt))
          (let ((l00 (list 'a 'b))
                (l01 (list 'c))
                (l02 (list 'ab 'ac))
                (l03 (list 'bc))
                (l10 (list 'c 'ab 'bc))
                (l11 (list 'ac))
                (l12 (list 'd 'bd 'cd))
                (l13 (list 'e 'be))
                (l20 (list 'd))
                (l21 (list 'e 'bd 'be 'de))
                (l23 (list 'bde))
                (l30 (list 'ac 'bd 'cd))
                (l31 (list 'abc))
                (l32 (list 'bde))
                (l33 (list 'bcd))
                (rslt +empty-list+))
            (declare (list l00 l01 l10 rslt))
            (progn
              (if (contains l00 gnrt)  (push (list '(0 0)) rslt))
              (if (contains l01 gnrt)  (push (list '(0 1)) rslt))
              (if (contains l02 gnrt)  (push (list '(0 2)) rslt))
              (if (contains l03 gnrt)  (push (list '(0 3)) rslt))
              (if (contains l10 gnrt)  (push (list '(1 0)) rslt))
              (if (contains l11 gnrt)  (push (list '(1 1)) rslt))
              (if (contains l12 gnrt)  (push (list '(1 2)) rslt))
              (if (contains l13 gnrt)  (push (list '(1 3)) rslt))
              (if (contains l20 gnrt)  (push (list '(2 0)) rslt))
              (if (contains l21 gnrt)  (push (list '(2 1)) rslt))
              (if (contains l23 gnrt)  (push (list '(2 3)) rslt))
              (if (contains l30 gnrt)  (push (list '(3 0)) rslt))
              (if (contains l31 gnrt)  (push (list '(3 1)) rslt))
              (if (contains l32 gnrt)  (push (list '(3 2)) rslt))
              (if (contains l33 gnrt)  (push (list '(3 3)) rslt))              
              (nreverse rslt)))))
;; #<Interpreted Function (unnamed) @ #x21304faa>

(setf ss4 (change-chcm-to-gflcc ss4 (dz2) ss-gflin4 'ss-gflin4))
;; [K31 Generalized-Filtered-Chain-Complex]

(gen-spsq-group ss4 (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((1 0) (0 1)),((1 1)),((1 1))]_{1}
;; Component Z

(gen-spsq-group ss4 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 1)),((1 1))]_{1}
;; Component Z
;; Component Z

(gen-spsq-group ss4 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 0) (0 1)),((1 1))]_{1}
;; Component Z

(gen-spsq-gnrts ss4 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
;; (
;; ----------------------------------------------------------------------{CMBN 1}
;; <1 * AB>
;; ------------------------------------------------------------------------------
;; 
;; ----------------------------------------------------------------------{CMBN 1}
;; <-1 * AC>
;; <1 * BC>
;; ------------------------------------------------------------------------------
;; )

(gen-spsq-gnrts ss4 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; (
;; ----------------------------------------------------------------------{CMBN 1}
;; <1 * AB>
;; ------------------------------------------------------------------------------
;; )

(gen-spsq-group ss4 (list '(2 1) '(1 2)) (list '(2 1) '(1 2)) (list '(3 3)) (list '(3 3)) 2)
;; Generalized spectral sequence S[((2 1) (1 2)),((2 1) (1 2)),((3 3)),((3 3))]_{2}
;; Component Z
;; Component Z
;; Component Z

(gen-spsq-gnrts ss4 (list '(2 1) '(1 2)) (list '(2 1) '(1 2)) (list '(3 3)) (list '(3 3)) 2)
;; (
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * ABC>
;; ------------------------------------------------------------------------------
;;  
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * BCD>
;; ------------------------------------------------------------------------------
;;  
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * BDE>
;; ------------------------------------------------------------------------------
;; )


;;; EXAMPLES OF CARTESIAN PRODUCTS OF THREE SIMPLICIAL SETS (TOWER OF FIBRATIONS WITH TRIVIAL TWISTING OPERATORS)


;; GENERALIZED FILTRATION FOR THE CARTESIAN PRODUCT


(setf crpr2-flin 
      #'(lambda (degr crpr)
          (declare
           (type fixnum degr)
           (type crpr crpr))
          (the fixnum
            (- degr (logcount (dgop1 crpr))))))
;; #<Interpreted Function (unnamed) @ #x20dd948a>


(setf crpr3-gflin 
      #'(lambda (degr gnrt)
          (declare 
           (type fixnum degr)
           (type gnrt gnrt))
          (let* ((dgop1 (dgop1 gnrt))
                 (gnrt1 (gnrt1 gnrt))
                 (flin1 (- degr (+ (logcount dgop1) (logcount (dgop1 gnrt1)))))
                 (flin2 (- degr (+ (logcount dgop1) (logcount (dgop2 gnrt1))))))
              (list (list (list flin1 flin2))))))
;; #<Interpreted Function (unnamed) @ #x20de254a>


;; GENERALIZED FILTRATION OF THE TENSOR PRODUCT (WHICH IS THE EFFECTIVE CHAIN COMPLEX OF THE CARTESIAN PRODUCT)



(setf tnpr2-flin 
      #'(lambda (degr tnpr)
          (declare 
           (ignore degr)
           (type tnpr tnpr))
          (the fixnum
            (degr1 tnpr))))
;; #<Interpreted Function (unnamed) @ #x20de85a2>


(setf tnpr3-gflin 
      #'(lambda (degr gnrt)
          (declare 
           (ignore degr)
           (type tnpr gnrt))
          (let ((flin1 (degr1 (gnrt1 gnrt)))
                (flin2 (degr2 (gnrt1 gnrt))))
            (list (list (list flin1 flin2))))))
;; #<Interpreted Function (unnamed) @ #x20deec42>


;;; EXAMPLE: S2 X S3 X S4

(setf s2 (sphere 2) s3 (sphere 3) s4 (sphere 4))
;; [K47 Simplicial-Set]

(setf x (crts-prdc (crts-prdc s2 s3) s4))
;; [K57 Simplicial-Set]

(setf xf (change-chcm-to-gflcc x (dz2) crpr3-gflin 'crpr3-gflin))
;; [K62 Generalized-Filtered-Chain-Complex]

(setf ex (rbcc (efhm xf)))
;; [K119 Chain-Complex]

(setf exf (change-chcm-to-gflcc ex (dz2) tnpr3-gflin 'tnpr3-gflin))
;; [K137 Generalized-Filtered-Chain-Complex]

(gfltrcm-efcc-fltr xf exf)
;; [K137 Generalized-Filtered-Chain-Complex]

(dotimes (n 6)
  (dotimes (p1 (1+ n))
    (dotimes (p2 (1+ (- n p1)))
      (if (e2-gspsq-gnrts xf p1 p2 n)
          (progn
            (print (list p1 p2))       
            (e2-gspsq-group xf p1 p2 n))))))
;; (0 0) Generalized spectral sequence S[((1 -2)),((1 -1)),((0 0)),((0 1) (1 0))]_{0}
;; Component Z
;; (2 0) Generalized spectral sequence S[((0 0) (1 -1) (3 -2)),((0 1) (1 0) (3 -1)),((0 1) (2 0)),((0 2) (2 1) (3 0))]_{2}
;; Component Z
;; (0 3) Generalized spectral sequence S[((1 1) (2 0)),((1 2) (2 1) (3 0)),((0 3) (1 2) (2 1) (3 0)),((0 4) (1 3) (2 2) (3 1) (4 0))]_{3}
;; Component Z
;; (0 0) Generalized spectral sequence S[((1 -2)),((1 -1)),((0 0)),((0 1) (1 0))]_{4}
;; Component Z
;; (2 3) Generalized spectral sequence S[((0 3) (1 2) (3 1) (4 0)),((0 4) (1 3) (3 2) (4 1) (5 0)),((0 4) (2 3) (3 2) (4 1) (5 0)),((0 5) (2 4) (3 3) (4 2) (5 1) (6 0))]_{5}
;; Component Z



;;; EXAMPLE: K(Z2,2) X K(Z2,3) X K(Z2,4)

(setf kz22 (k-z2 2) kz23 (k-z2 3) kz24 (k-z2 4))
;; [K175 Abelian-Simplicial-Group]

(setf y (crts-prdc (crts-prdc kz22 kz23) kz24))
;; [K192 Simplicial-Set]

(setf yf (change-chcm-to-gflcc y (dz2) crpr3-gflin 'crpr3-gflin))
;; [K197 Generalized-Filtered-Chain-Complex]

(setf ey (rbcc (efhm y)))
;; [K612 Chain-Complex]

(setf eyf (change-chcm-to-gflcc ey (dz2) tnpr3-gflin 'tnpr3-gflin))
;; [K633 Generalized-Filtered-Chain-Complex]

(gfltrcm-efcc-fltr yf eyf)
;; [K633 Generalized-Filtered-Chain-Complex]

;; BY DEFAULT THE COMPUTATIONS ARE DONE BY USING EFFECTIVE HOMOLOGY
(dotimes (n 6)
  (dotimes (p1 (1+ n))
    (dotimes (p2 (1+ (- n p1)))
      (if (e2-gspsq-gnrts yf p1 p2 n)
          (progn
            (print (list p1 p2))       
            (e2-gspsq-group yf p1 p2 n))))))

;; (0 0) Generalized spectral sequence S[((1 -2)),((1 -1)),((0 0)),((0 1) (1 0))]_{0}
;; Component Z
;; (2 0) Generalized spectral sequence S[((0 0) (1 -1) (3 -2)),((0 1) (1 0) (3 -1)),((0 1) (2 0)),((0 2) (2 1) (3 0))]_{2}
;; Component Z/2Z
;; (0 3) Generalized spectral sequence S[((1 1) (2 0)),((1 2) (2 1) (3 0)),((0 3) (1 2) (2 1) (3 0)),((0 4) (1 3) (2 2) (3 1) (4 0))]_{3}
;; Component Z/2Z
;; (0 0) Generalized spectral sequence S[((1 -2)),((1 -1)),((0 0)),((0 1) (1 0))]_{4}
;; Component Z/2Z
;; (4 0) Generalized spectral sequence S[((0 2) (1 1) (2 0) (3 -1) (5 -2)),((0 3) (1 2) (2 1) (3 0) (5 -1)),((0 3) (1 2) (2 1) (4 0)),((0 4) (1 3) (2 2) (4 1) (5 0))]_{4}
;; Component Z/4Z
;; (0 5) Generalized spectral sequence S[((1 3) (2 2) (3 1) (4 0)),((1 4) (2 3) (3 2) (4 1) (5 0)),((0 5) (1 4) (2 3) (3 2) (4 1) (5 0)),((0 6) (1 5) (2 4) (3 3) (4 2) (5 1) (6 0))]_{5}
;; Component Z/2Z
;; (2 3) Generalized spectral sequence S[((0 3) (1 2) (3 1) (4 0)),((0 4) (1 3) (3 2) (4 1) (5 0)),((0 4) (2 3) (3 2) (4 1) (5 0)),((0 5) (2 4) (3 3) (4 2) (5 1) (6 0))]_{5}
;; Component Z/2Z
;; (5 0) Generalized spectral sequence S[((0 3) (1 2) (2 1) (3 0) (4 -1) (6 -2)),((0 4) (1 3) (2 2) (3 1) (4 0) (6 -1)),((0 4) (1 3) (2 2) (3 1) (5 0)),((0 5) (1 4) (2 3) (3 2) (5 1) (6 0))]_{5}
;; Component Z/2Z



;; COMPUTATIONS WITHOUT EFFECTIVE HOMOLOGY ARE ALSO POSSIBLE (THE RESULTS ARE THE SAME BUT ONLY POSSIBLE UNTIL DEGREE 3)


(dotimes (n 4)
  (dotimes (p1 (1+ n))
    (dotimes (p2 (1+ (- n p1)))
      (let* ((point (list p1 p2))
             (p (t-downset-list point))
             (q (t-downset-list (2-points-add point (list 1 -1))))
             (z (t-downset-list (2-points-add point (list 1 -2))))
             (b (t-downset-list (2-points-add point (list 0 1)))))
        (declare (type list point p q z b))
        (if (gen-eff-spsq-gnrts yf z q p b n)
            (progn
              (print point)
              (gen-eff-spsq-group yf z q p b n)))))))

;; (0 0) Generalized spectral sequence S[((1 -2)),((1 -1)),((0 0)),((0 1) (1 0))]_{0}
;; Component Z
;; (2 0) Generalized spectral sequence S[((0 0) (1 -1) (3 -2)),((0 1) (1 0) (3 -1)),((0 1) (2 0)),((0 2) (2 1) (3 0))]_{2}
;; Component Z/2Z
;; (0 3) Generalized spectral sequence S[((1 1) (2 0)),((1 2) (2 1) (3 0)),((0 3) (1 2) (2 1) (3 0)),((0 4) (1 3) (2 2) (3 1) (4 0))]_{3}
;; Component Z/2Z