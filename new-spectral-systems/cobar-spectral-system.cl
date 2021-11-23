;;;  COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM
;;;  COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM
;;;  COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM COBAR-SPECTRAL-SYSTEM


(IN-PACKAGE #:cat)

(PROVIDE "cobar-spectral-system")


(DEFVAR TNPR-COBAR-GFLIN)

(SETF TNPR-COBAR-GFLIN 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((degr1 (degr1 gnrt))
            (allp (gnrt1 gnrt))
            (degr2 (degr2 gnrt))
            (tnpr2 (gnrt2 gnrt)))
        (declare (ignore degr1))
        (let ((flin1 (funcall tnpr-flin degr2 tnpr2))
              (flin2 (with-allp (l) allp
                       (the fixnum
                         (- (length l))))))
          (list (list (list flin1 flin2)))))))

(DEFUN CHECK-TNPR-COBAR-GFLIN (fltrcm n1 &rest rest)
  (let ((gf (change-chcm-to-gflcc fltrcm (dz2) tnpr-cobar-gflin '(tnpr-cobar-gflin))))
    (let ((n2 (if (> (length rest) 0) (first rest) (1+ n1))))
      (dotimes (n0 (- n2 n1))
        (let ((n (+ n0 n1)))
          (dotimes (i (length (basis gf n)))
            (if (not (gfltrd-mrph-order-p (dffr gf) n (nth i (basis gf n)) '(0 0)))
                (return-from check-tnpr-cobar-gflin nil)))))
      t)))


(DEFVAR TNPR-COBAR-GFLIN2)

(SETF TNPR-COBAR-GFLIN2 
  #'(lambda (degr gnrt)
      (declare 
       (ignore degr)
       (type tnpr gnrt))
      (let ((degr1 (degr1 gnrt))
            (allp (gnrt1 gnrt))
            (degr2 (degr2 gnrt))
            (tnpr2 (gnrt2 gnrt)))
        (declare (ignore tnpr2 degr1))
        (let ((flin1 degr2)
              (flin2 (with-allp (l) allp
                       (the fixnum
                         (- (length l))))))
          (list (list (list flin1 flin2)))))))

(DEFUN CHECK-TNPR-COBAR-GFLIN2 (fltrcm n1 &rest rest)
  (let ((gf (change-chcm-to-gflcc fltrcm (dz2) tnpr-cobar-gflin2 '(tnpr-cobar-gflin2))))
    (let ((n2 (if (> (length rest) 0) (first rest) (1+ n1))))
      (dotimes (n0 (- n2 n1))
        (let ((n (+ n0 n1)))
          (dotimes (i (length (basis gf n)))
            (if (not (gfltrd-mrph-order-p (dffr gf) n (nth i (basis gf n)) '(0 0)))
                (return-from check-tnpr-cobar-gflin2 nil)))))
      t)))


(DEFUN COBAR-SPECTRAL-SYSTEM (fibration)
  (let* ((k (eilenberg-moore-bicomplex fibration))
         (ecc (rbcc (efhm k))))
    (declare 
     (type chain-complex k)
     (type chain-complex ecc))
    (progn
      (setf k (change-chcm-to-gflcc k (dz2) tnpr-cobar-gflin '(tnpr-cobar-gflin)))
      (setf ecc (change-chcm-to-gflcc ecc (dz2) tnpr-cobar-gflin '(tnpr-cobar-gflin)))
      (the SPECTRAL-SYSTEM
        (build-spectral-system ecc `(cobar-spectral-system ,fibration))))))  






