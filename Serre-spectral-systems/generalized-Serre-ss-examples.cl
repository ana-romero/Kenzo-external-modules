(in-package :user)
(load "cat-init")
(load-cfiles)
(load "filtered-complexes.cl")
(load "spectral-sequences.lsp")
(load "generalized-spectral-sequences.cl")
(load "generalized-serre-ss.cl")


;;; POSTNIKOV TOWER OF THE SPHERE S^3 WITH TWO FIBRATIONS
(CAT-INIT)

(progn
  (setf B (sphere 3))
  (setf k1 (chml-clss B 3))
  (setf t1 (z-whitehead B k1))
  (setf N (fibration-total t1))
  (setf k0 (chml-clss N 4))
  (setf t0 (z2-whitehead N k0))
  (setf E (fibration-total t0)))
;; [K298 Simplicial-Set]


;;; EFFECTIVE HOMOLOGY OF E

(efhm E)
;; [K608 Homotopy-Equivalence K298 <= K598 => K594]

;;; DEFINITION OF THE FILTRATIONS ON THE SPACES

(setf Ef (change-chcm-to-gflcc E (dz2) crpr2-gflin 'crpr2-gflin))
;; [K611 Generalized-Filtered-Chain-Complex]

(setf Df (change-chcm-to-gflcc (rbcc (efhm e)) (dz2) tnpr2-gflin 'tnpr2-gflin))
;; [K613 Generalized-Filtered-Chain-Complex]

(gfltrcm-efcc-fltr Ef Df)

(setf E ef)

;;; COMPUTATION OF THE LEVEL 2 OF THE SERRE SPECTRAL SYSTEM

(e2-gspsq-group E '(0 0) 5)
;; Component Z/2Z

(e2-gspsq-group E '(2 0) 5)
;; Component Z/2Z

(e2-gspsq-group E '(2 3) 5)
;; Component Z


;;; LEXICOGRAPHICAL CONNECTIONS

(lexcon-gspsq-group E '(0 0) '(1 1) 5)
;; Component Z/2Z

(lexcon-gspsq-group E '(2 0) '(1 1) 5)
;; NIL

(lexcon-gspsq-group E '(2 3) '(1 1) 5)
;; NIL



;;; FINAL GROUP OF THE SPECTRAL SYSTEM FOR DEGREE 5

(final-gspsq-group E 5)
;; Component Z/2Z



;;; POSTNIKOV TOWER OF THE SPHERE S^3 WITH THREE FIBRATIONS



(cat-init)
(progn
  (setf s3 (sphere 3))
  (setf k3 (chml-clss s3 3))
  (setf F3 (z-whitehead s3 k3))
  (setf X4 (fibration-total F3))
  (setf k4 (chml-clss x4 4))
  (setf F4 (z2-whitehead X4 k4))
  (setf X5 (fibration-total F4))
  (setf k5 (chml-clss x5 5))
  (setf F5 (z2-whitehead X5 k5))
  (setf X6 (fibration-total F5)))
;; [K630 Simplicial-Set]

;;; GENERALIZED FILTRATIONS

(setf X (change-chcm-to-gflcc X6 (downsets (zn 3)) crpr3-gflin 'crpr3-gflin))
;; [K637 Generalized-Filtered-Chain-Complex]

(setf XEf (change-chcm-to-gflcc (rbcc (efhm X)) (downsets (zn 3)) tnpr3-gflin 'tnpr3-gflin))
;; [K845 Generalized-Filtered-Chain-Complex]

(gfltrcm-efcc-fltr X XEf)
;; [K845 Generalized-Filtered-Chain-Complex]

;;; LEVEL 2 OF THE SPECTRAL SYSTEM

(e2-gspsq-group X '(0 6 0) 6)
;; Component Z

;;; LEXICOGRAPHICAL CONNECTION

(lexcon-gspsq-group X '(0 6 0) '(1 1 1) 6)
;; Component Z/3Z

;;; FINAL GROUP OF THE SPECTRAL SYSTEM FOR DEGREE 6

(final-gspsq-group X 6)
;; Component Z/12Z


;;; EFFECTIVE EXAMPLE

(cat-init)
(progn
  (setf kz22 (k-z2 2))
  (setf k2 (chml-clss kz22 4))
  (setf F2 (z2-whitehead kz22 k2))
  (setf X3 (fibration-total f2))
  (setf k3 (chml-clss x3 5))
  (setf F3 (z2-whitehead X3 k3))
  (setf X4 (fibration-total F3))
  (setf k4 (chml-clss x4 6))
  (setf F4 (z2-whitehead X4 k4))
  (setf X5 (fibration-total F4)))
;; [K612 Kan-Simplicial-Set]

(setf ex5 (rbcc (efhm x5)))
;; [K808 Chain-Complex]

(setf effX5 (change-chcm-to-gflcc ex5 (dzn 3) tnpr3-gflin 'tnpr3-gflin))
;; [K825 Generalized-Filtered-Chain-Complex]

(setf x5f (change-chcm-to-gflcc x5 (dzn 3) crpr3-gflin 'crpr3-gflin))
;; [K827 Generalized-Filtered-Chain-Complex]

(setf x5 x5f)
;;[K827 Generalized-Filtered-Chain-Complex]

(e1-eff-gspsq-group X5 '(0 0 3) 3)
;; Component Z
;; Component Z
;; Component Z
;; Component Z
;; Component Z
;; Component Z
;; Component Z
;; Component Z

(e1-eff-gspsq-group effX5 '(0 0 3) 3)
;; Component Z