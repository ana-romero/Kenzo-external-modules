(in-package :user)
(load "cat-init")
(load-cfiles)
(load "filtered-complexes.cl")
(load "spectral-sequences.cl")
(load "generalized-spectral-sequences.cl")
(load "chain-complex-vf")
(load "vf-order")
(load "generalized-filtered-vf")


;; WE CONSTRUCT THE SAME SIMPLICIAL SET ss1 AS IN THE FILE generalized-spectral-sequences-examples.cl

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


(setf mtrx-degr 1)
;; 1

(setf vf-list (gfltrcm-vf-list ss1 1))
;; ((((1 0)) NIL) (((0 1)) ((0 0))) (((0 1) (1 0)) ((0 0))) (((0 0)) ((0 0) (2 1))))
;; The function gfltrcm-vf-list returns a list of pairs. The first element is the generalized
;; filtration index (a list of elements of the poset) and the second one is a list of pairs defining
;; the vector field.
;; For example, the first pair is (((1 0)) NIL) and this means that for filtration index ((1 0)) there are not vectors.
;; The second pair is (((0 1)) ((0 0))) and this means that for filtration index ((0 1)) there is a vector (0 0)
;; which corresponds to tau=element 0 of filtration index ((0 1)) and sigma= element 0 of filtration index ((0 1)).


;;; LIST OF CRITICAL CELLS IN DEGREE 1
(gfltrcm-vf-crtc-list ss1 mtrx-degr vf-list 1)
;;(CD DE BC)


;;; LIST OF CRITICAL CELLS IN DEGREE 0
(gfltrcm-vf-crtc-list ss1 mtrx-degr vf-list 0)
(B)

(setf crtc-chcm (filtered-vf-crtc-chcm ss1 mtrx-degr vf-list))
;; [K11 Generalized-Filtered-Chain-Complex]


;; COMPUTATION OF THE GENERALIZED SPECTRAL SEQUENCE OF THE CRITICAL COMPLEX
;; (WE OBTAIN THE SAME RESULTS AS FOR ss1, AS EXPECTED)
(gen-spsq-group crtc-chcm '(0 0) '(0 0) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(0 0),(0 0),(1 1),(1 1)]_{1}
;; Component Z
;; Component Z

(gen-spsq-group crtc-chcm '(0 0) '(0 1) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(0 0),(0 1),(1 1),(1 1)]_{1}
;; Component Z

(gen-spsq-group crtc-chcm '(1 0) '(1 0) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(1 0),(1 0),(1 1),(1 1)]_{1}
;; Component Z

(gen-spsq-group crtc-chcm '(1 -1) '(1 -1) '(1 1) '(1 1) 1)
;; Generalized spectral sequence S[(1 -1),(1 -1),(1 1),(1 1)]_{1}
;; Component Z
;; Component Z
;; Component Z

(gen-spsq-group crtc-chcm '(0 0) '(0 0) '(1 1) '(1 1) 2)
;; Generalized spectral sequence S[(0 0),(0 0),(1 1),(1 1)]_{2}
;; NIL


;; WE CAN ALSO CONSTRUCT DIRECTLY THE REDUCTION ss1 ==> crtc-chcm
(setf r1 (gfltr-vf-1degr-rdct ss1 1))
[K21 Reduction K7 => K11]
;; We observe the effective chain complex in the reduction [K11] which is crtc-chcm


;; ANOTHER EXAMPLE A BIT MORE COMPLICATED IS ss4 (THE ONE IN THE FILE generalized-spectral-sequences-examples)
;; BECAUSE THERE ARE ALSO 2-SIMPLICES. THE FILTRATION IS DEFINED OVER DZ2


(setf ss4
  (build-finite-ss
   '(a b c d e
       1 ab (b a) ac (c a) bc (c b) bd (d b) cd (d c) be (e b) de (e d)
       2 abc (bc ac ab) bcd (cd bd bc) bde (de be bd)
        )))
;; [K22 Simplicial-Set]

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
;; [K27 Generalized-Filtered-Chain-Complex]


;; IN THIS CASE WE CAN CONSIDER TWO DIFFERENT VECTOR FIELDS, ONE FOR THE MATRIX OF DIMENSION 1
;; AND ANOTHER ONE FOR THE MATRIX OF DIMENSION 2. WE CAN DO IT SEPARATELY AS BEFORE OR DIRECTLY
;; WITH THE FOLLOWING FUNCTION

(setf r4 (gfltr-vf-rdct ss4 2))
;; [K54 Reduction K27 => K43]

;; IN THIS CASE THE REDUCTION IS BUILT AS THE COMPOSITION OF TWO REDUCTIONS, WITH INTERMEDIATE
;; COMPLEX K36
(orgn r4)
;; (CMPS [K48 Reduction K36 => K43] [K41 Reduction K27 => K36])


(setf crtc4 (bcc r4))
;; [K43 Generalized-Filtered-Chain-Complex]

;; THE RESULTS OF THE GENERALIZED SPECTRAL SEQUENCE ARE THE SAME AS FOR ss4
;; (ONLY THE ORDER OF SOME GENERATORS HAS CHANGED)
(gen-spsq-group crtc4 (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((1 0) (0 1)),((1 1)),((1 1))]_{1}
;; Component Z

(gen-spsq-group crtc4 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 1)),((1 1))]_{1}
;; Component Z
;; Component Z

(gen-spsq-group crtc4 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; Generalized spectral sequence S[((0 0)),((0 0)),((1 0) (0 1)),((1 1))]_{1}
;; Component Z

(gen-spsq-gnrts crtc4 (list '(0 0)) (list '(0 0)) (list '(1 1)) (list '(1 1)) 1)
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

(gen-spsq-gnrts crtc4 (list '(0 0)) (list '(0 0)) (list '(1 0) '(0 1)) (list '(1 1)) 1)
;; (
;; ----------------------------------------------------------------------{CMBN 1}
;; <1 * AB>
;; ------------------------------------------------------------------------------
;; )

(gen-spsq-group crtc4 (list '(2 1) '(1 2)) (list '(2 1) '(1 2)) (list '(3 3)) (list '(3 3)) 2)
;; Generalized spectral sequence S[((2 1) (1 2)),((2 1) (1 2)),((3 3)),((3 3))]_{2}
;; Component Z
;; Component Z
;; Component Z

(gen-spsq-gnrts crtc4 (list '(2 1) '(1 2)) (list '(2 1) '(1 2)) (list '(3 3)) (list '(3 3)) 2)
;; (
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * BDE>
;; ------------------------------------------------------------------------------
;;  
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * BCD>
;; ------------------------------------------------------------------------------
;;  
;; ----------------------------------------------------------------------{CMBN 2}
;; <1 * ABC>
;; ------------------------------------------------------------------------------
;; )




