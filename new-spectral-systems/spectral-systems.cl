;;;  SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS
;;;  SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS
;;;  SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS SPECTRAL-SYSTEMS


(IN-PACKAGE #:cat)

(PROVIDE "spectral-systems")

(DEFCLASS SPECTRAL-SYSTEM ()
  (;; GeneralizedFiLTeRedchainCoMplex
   (gfltrcm :type generalized-filtered-chain-complex :initarg :gfltrcm :reader gfltrcm)
   ;; ReSuLTS 
   (group-rslts :type list :initarg :group-rslts :reader group-rslts)
   (dffr-rslts :type list :initarg :dffr-rslts :reader dffr-rslts)
   ;; IDentification NuMber      
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)      
   ;; ORiGiN      
   (orgn :type list :initarg :orgn :reader orgn)))

(DEFVAR *ssy-list*
    "The variable *SSY-LIST* is bound to a list of user created spectral systems")
(SETF *ssy-list* +empty-list+)
(PUSHNEW '*ssy-list* *list-list*)


#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) nil))


(DEFMETHOD PRINT-OBJECT ((ss SPECTRAL-SYSTEM) stream)
  (the spectral-system
    (progn
      (format stream "[K~D Spectral-System]" (idnm ss))
      ss)))

#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) t))


(DEFUN SSY (idnm)
  (declare (type fixnum idnm))
  (the (or spectral-system null)
    (find idnm *ssy-list* :key #'idnm)))


;;; Function to build a spectral systems from a generalized filtered chain complex and an origin
(DEFUN BUILD-SPECTRAL-SYSTEM (gfltrcm orgn)
  (declare
   (type generalized-filtered-chain-complex gfltrcm)
   (type list orgn))
  (the spectral-system
    (progn
      (let ((already (find orgn *ssy-list* :test #'equal :key #'orgn)))
        (declare (type (or spectral-system null) already))
        (when already
          (return-from build-spectral-system already)))
      (let ((ss (make-instance 'spectral-system
                  :gfltrcm gfltrcm
                  :orgn orgn
                  :group-rslts +empty-list+
                  :dffr-rslts +empty-list+)))
        (declare (type spectral-system ss))
        (push ss *ss-list*)
        ss))))


;; Method that returns the representation basis-divisors of the group S[z,p,s,b]_n of the spectral system
(DEFMETHOD SPECTRAL-SYSTEM-BASIS-DVS ((ss SPECTRAL-SYSTEM) z s p b degr)
  (declare 
   (type spectral-system ss)
   (type fixnum degr)
   (type gnrt z s p b))
  (let* ((gfltrcm (gfltrcm ss))
         (rslts (group-rslts ss))
         (pos (position (list z s p b degr nil) rslts :test #'(lambda (l1 l2)
                                                         (and
                                                          (eq (first l1) (first l2))
                                                          (eq (second l1) (second l2))
                                                          (eq (third l1) (third l2))
                                                          (eq (fourth l1) (third l2))
                                                          (eq (fifth l1) (third l2)))))))
    (the list
      (if pos (sixth (nth pos rslts))
        (let ((basis-dvs (GEN-SPSQ-BASIS-DVS gfltrcm z s p b degr)))
          (declare (type list basis-dvs))
          (progn
            (push (list z s p b degr basis-dvs) rslts)
            (setf (slot-value ss 'group-rslts) rslts)
            basis-dvs))))))


;; Method that returns the list of abelian invariants of the group S[z,p,s,b]_n of the spectral system
(DEFMETHOD SPECTRAL-SYSTEM-GROUP ((ss SPECTRAL-SYSTEM) z s p b degr)
  (declare 
   (type spectral-system ss)
   (type fixnum degr)
   (type gnrt z s p b))  
  (let* ((basis-dvs (spectral-system-basis-dvs ss z s p b degr))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (mapcan #'(lambda (i)
                (if (eq 1 i) nil
                  (list i)))
      divs)))


;; Method that prints the spectral system group
(DEFMETHOD PRINT-SPECTRAL-SYSTEM-GROUP ((ss SPECTRAL-SYSTEM) z s p b degr)
  (declare 
   (type spectral-system ss)
   (type fixnum degr)
   (type gnrt z s p b))  
  (let* ((gfltrcm (gfltrcm ss)))
    (gen-spsq-group gfltrcm z s p b degr)))


;; Method that determines de differential map of a spectral system applied on an element given by its coordinates
;; in int-list
(DEFMETHOD SPECTRAL-SYSTEM-DIFFERENTIAL ((ss SPECTRAL-SYSTEM) z1 s1 p1 b1 z2 s2 p2 b2 degr int-list)
  (declare
   (type spectral-system ss)
   (type fixnum degr)
   (type gnrt z1 s1 p1 b1 z2 s2 p2 b2)
   (type list int-list))
  (let* ((gfltrcm (gfltrcm ss)))
    (gen-spsq-dffr gfltrcm z1 s1 p1 b1 z2 s2 p2 b2 degr int-list)))


