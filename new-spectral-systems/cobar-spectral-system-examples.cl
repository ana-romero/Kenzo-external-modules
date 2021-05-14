;;; 
(load "Kenzo-external-modules\\Eilenberg-Moore-spectral-sequence\\central-extensions")
(load "Kenzo-external-modules\\Eilenberg-Moore-spectral-sequence\\eilenberg-moore-spectral-sequence-auxfunctions")
(load "Kenzo-external-modules\\Eilenberg-Moore-spectral-sequence\\eilenberg-moore-spectral-sequence")
(load "Kenzo-external-modules\\generalized-spectral-sequences\\generalized-spectral-sequences.cl")
(load "Kenzo-external-modules\\Serre-spectral-systems\\generalized-Serre-ss.cl")
(load "Kenzo-external-modules\\new-spectral-systems\\spectral-systems")
(load "Kenzo-external-modules\\new-spectral-systems\\cobar-spectral-system")


(IN-PACKAGE #:cat)


;; Example 1: Cobar spectral system of the fibration K(Z/4,1) -> K(Z/2,1) -> K(Z/2, 2) corresponding to the
;; extension Z/2 -> Z/4 -> Z/2

(cat-init) 
(setf A (cyclicGroup 2))
(setf B (cyclicGroup 2))
(setf KA1 (k-g-1 A))
(setf KB1 (k-g-1 B))                    

(setf cocycle #'(lambda (z1 z2)
                  (if (and (= z1 1) (= z2 1))
                      1 0)))

(setf E (gr-cntr-extn a b cocycle))
(change-class E 'AB-group)

(setf KE1 (k-g-1 E))
(setf efhm2 (central-extension-efhm A B cocycle))
(setf (slot-value kE1 'efhm) efhm2)
(setf ka2 (k-g a 2))
(setf tau1 (univ-fbrt-tw A 2))
(setf tw
  (the fibration
    (build-smmr :sorc ka2 :trgt ke1 :degr -1
                :sintr 
                (sintr (cmps (cocycle-fibr-iso2 b a cocycle)
                             (cmps  (twop-incl (cocycle-fibration B A cocycle)) tau1)))
                :orgn `(fibration33 ,tau1))))


(setf ss1 (cobar-spectral-system tw))

(PRINT-SPECTRAL-SYSTEM-GROUP ss1 (list '(0 0)) (list '(0 0)) (list '(1 0)) (list '(1 0)) 4)

(setf r 2)
(dotimes (n 5)
  (dotimes (p (1+ n))
    (PRINT-SPECTRAL-SYSTEM-GROUP ss1 (list (list (- p r) 0))  (list (list (- p 1) 0))  (list  (list p 0))  (list (list (- (+ p r) 1) 0)) n)))



;; Example 2: 

(cat-init)
(setf A (z-group))
(setf B (gr-crts-prdc (z-group) (z-group)))
(setf cocycle #'(lambda (crpr1 crpr2)
                  (with-grcrpr (x1 y1) crpr1
                    (with-grcrpr (x2 y2) crpr2
                      (div (- (* x1 y2) (* y1 x2)) 2)))))
(setf E (gr-cntr-extn A B cocycle))

(setf KA1 (k-g-1 A))
(efhm KA1)
(setf KB1 (k-g-1 B))
(setf efhm (gr-crts-prdc-efhm A A))
(setf (slot-value kB1 'efhm) efhm)

(setf KE1 (k-g-1 E))
(setf efhm2 (central-extension-efhm A B cocycle))
(setf (slot-value kE1 'efhm) efhm2)

(setf ka2 (k-g a 2))

(setf tau1 (univ-fbrt-tw A 2))

(setf tw
  (the fibration
    (build-smmr :sorc ka2 :trgt ke1 :degr -1
                :sintr 
                (sintr (cmps (cocycle-fibr-iso2 b a cocycle)
                             (cmps  (twop-incl (cocycle-fibration B A cocycle)) tau1)))
                :orgn `(fibration33 ,tau1))))

(setf ss2 (cobar-spectral-system tw))

(PRINT-SPECTRAL-SYSTEM-GROUP ss2 (list '(0 0)) (list '(0 0)) (list '(1 0)) (list '(1 0)) 4)

(setf r 2)
(dotimes (n 5)
  (dotimes (p (1+ n))
    (PRINT-SPECTRAL-SYSTEM-GROUP ss2 (list (list (- p r) 0))  (list (list (- p 1) 0))  (list  (list p 0))  (list (list (- (+ p r) 1) 0)) n)))


