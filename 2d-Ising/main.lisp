;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;Global Variables and Structures;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *N* 10)

(defstruct spin
  value-x
  value-y
  value-z)

(defparameter *3d-spin-lattice*
  (make-array (list *N* *N* *N*)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; MACROS ;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro do-lattice ((spin x y z lattice) &body body)
  "Iterates over a 3D lattice.
   
   Arguments:
     SPIN: Symbol to bind to the current spin struct.
     X, Y, Z: Symbols to bind to the current coordinates.
     LATTICE: The 3D array to iterate over.
     BODY: The code to execute for each spin.
   
   Returns:
     NIL (implicitly via dotimes)."
  (let ((l-var (gensym))
        (dim-x (gensym)) (dim-y (gensym)) (dim-z (gensym)))
    `(let* ((,l-var ,lattice)
            (,dim-x (array-dimension ,l-var 0))
            (,dim-y (array-dimension ,l-var 1))
            (,dim-z (array-dimension ,l-var 2)))
       (dotimes (,x ,dim-x)
         (dotimes (,y ,dim-y)
           (dotimes (,z ,dim-z)
             (let ((,spin (aref ,l-var ,x ,y ,z)))
               ,@body)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Periodic Boundary Helper Function & Neighbors;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pbc (index max-index)
  "Wraps an index to fit within [0, max-index).
   Used for periodic boundary conditions."
  (mod index max-index))

(defun spin-at (lattice i j k)
  "Safe accessor for the lattice.
   
   Arguments:
     LATTICE: The spin array.
     I, J, K: Integer coordinates (can be negative or > N).
   
   Returns:
     The spin struct at (i,j,k) wrapped periodically."
  (let ((N (array-dimension lattice 0)))
    (aref lattice 
          (pbc i N) 
          (pbc j N) 
          (pbc k N))))

(defun get-neighbors (lattice i j k)
  "Arguments:
     LATTICE: The spin array.
     I, J, K: Center coordinates.
   
   Returns:
     List of 6 neighboring spin structs (+x, -x, +y, -y, +z, -z)."
  (list (spin-at lattice (1+ i) j k)
        (spin-at lattice (1- i) j k)
        (spin-at lattice i (1+ j) k)
        (spin-at lattice i (1- j) k)
        (spin-at lattice i j (1+ k))
        (spin-at lattice i j (1- k))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;Spin Observables;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun spin-magnitude (s)
  "Calculates the Euclidean norm of the spin vector."
  (sqrt (+ (expt (spin-value-x s) 2)
           (expt (spin-value-y s) 2)
           (expt (spin-value-z s) 2))))

(defun spin-angles (s)
  "Converts Cartesian spin components to Spherical coordinates.
   
   Arguments:
     S: Spin struct.
   
   Returns (2 Values):
     THETA: Polar angle (radians, 0 to pi).
     PHI: Azimuthal angle (radians, -pi to pi)."
  (let ((mag (spin-magnitude s)))
    (if (< mag 1.0e-6)
        (values 0.0 0.0)
        (values 
         (acos (/ (spin-value-z s) mag))
         (atan (spin-value-y s) (spin-value-x s))))))

(defun calculate-magnetization (lattice)
  "Calculates the system's total magnetization order parameter.
   
   Arguments:
     LATTICE: The spin array.
   
   Returns (4 Values):
     1. AVG-SPIN: A spin struct representing the average vector.
     2. MAG: The magnitude |M| (0.0 to 1.0).
     3. THETA: The polar angle of the global moment.
     4. PHI: The azimuthal angle of the global moment."
  (let ((total-spin (make-spin :value-x 0.0 :value-y 0.0 :value-z 0.0))
        (N (array-total-size lattice)))
    
    (do-lattice (s i j k lattice)
      (incf (spin-value-x total-spin) (spin-value-x s))
      (incf (spin-value-y total-spin) (spin-value-y s))
      (incf (spin-value-z total-spin) (spin-value-z s)))
    
    (setf (spin-value-x total-spin) (/ (spin-value-x total-spin) N))
    (setf (spin-value-y total-spin) (/ (spin-value-y total-spin) N))
    (setf (spin-value-z total-spin) (/ (spin-value-z total-spin) N))
    
    (let ((mag (spin-magnitude total-spin)))
      (multiple-value-bind (theta phi) (spin-angles total-spin)
        (values total-spin mag theta phi)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;Lattice Initialization;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun random-unit-vector ()
  "Generates a random unit vector uniformly distributed on a sphere.
   Returns: Values X, Y, Z."
  (let* ((phi (random (* 2.0 pi)))
         (z   (- (random 2.0) 1.0))
         (r   (sqrt (- 1.0 (* z z)))))
    (values (* r (cos phi))
            (* r (sin phi))
            z)))

(defun normalize-spin! (s)
  "Normalizes a spin struct in-place to have magnitude 1.0."
  (let ((mag (spin-magnitude s)))
    (when (> mag 0.0)
      (setf (spin-value-x s) (/ (spin-value-x s) mag))
      (setf (spin-value-y s) (/ (spin-value-y s) mag))
      (setf (spin-value-z s) (/ (spin-value-z s) mag))))
  s)

(defun initialize-spins (lattice &optional (randomly t))
  "Resets the lattice spins.
   
   Arguments:
     LATTICE: The array to modify.
     RANDOMLY: If T, directions are random. If NIL, all point Up (+Z)."
  (do-lattice (s i j k lattice)
    (declare (ignore s))
    (let ((new-spin 
            (if randomly
                (multiple-value-bind (s-x s-y s-z) (random-unit-vector)
                  (make-spin :value-x s-x :value-y s-y :value-z s-z))
                (make-spin :value-x 0.0 :value-y 0.0 :value-z 1.0))))
      (setf (aref lattice i j k) new-spin)
      (normalize-spin! new-spin))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;Calculation of Hamiltonian;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun spin-dot-product (s1 s2)
  "Calculates the scalar dot product of two spin vectors."
  (+ (* (spin-value-x s1) (spin-value-x s2))
     (* (spin-value-y s1) (spin-value-y s2))
     (* (spin-value-z s1) (spin-value-z s2))))

(defun calculate-hamiltonian (lattice &optional (Js '(1.0 1.0 1.0 1.0 1.0 1.0)))
  "Calculates the total energy of the system.
   
   Arguments:
     LATTICE: The spin array.
     JS: List of 6 coupling constants (+x -x +y -y +z -z).
   
   Returns:
     Total Energy / 2.0 (Float)."
  (let ((total-energy 0.0))
    (do-lattice (s1 i j k lattice)
      (let ((neighbors (get-neighbors lattice i j k)))
        (mapc (lambda (s2 J) 
                (decf total-energy (* J (spin-dot-product s1 s2))))
              neighbors Js)))
    (/ total-energy 2.0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;Metropolis Algorithm & Simulation;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun metropolis-sweep (lattice Temp &optional (Js '(1.0 1.0 1.0 1.0 1.0 1.0)))
  "Performs one Monte Carlo sweep (attempting to flip every spin once).
   
   Arguments:
     LATTICE: The spin array.
     TEMP: Temperature (kT).
     JS: Coupling constants.
   
   Returns:
     Integer: Number of accepted flips."
  (let ((beta (if (> Temp 0.0) (/ 1.0 Temp) 10000000000.0))
        (accepted-flips 0))
    
    (do-lattice (s i j k lattice)
      (let ((neighbors (get-neighbors lattice i j k)))
        (multiple-value-bind (nx ny nz) (random-unit-vector)
          
          (let ((delta-E 0.0)
                (dx (- nx (spin-value-x s)))
                (dy (- ny (spin-value-y s)))
                (dz (- nz (spin-value-z s))))
            
            (mapc (lambda (n J)
                    (let ((dot-product (+ (* (spin-value-x n) dx)
                                          (* (spin-value-y n) dy)
                                          (* (spin-value-z n) dz))))
                      (decf delta-E (* J dot-product))))
                  neighbors Js)
            
            (when (or (< delta-E 0.0)
                      (< (random 1.0) (exp (* (- delta-E) beta))))
              (setf (spin-value-x s) nx)
              (setf (spin-value-y s) ny)
              (setf (spin-value-z s) nz)
              (incf accepted-flips))))))
    accepted-flips))

(defun run-simulation (steps temp &optional (Js '(1.0 1.0 1.0 1.0 1.0 1.0)) (cold-steps 0) (record-interval 1) (print-interval 100) (print-summary t))
  "Runs the Metropolis simulation.
   
   Arguments:
     STEPS: Number of Monte Carlo steps to record.
     TEMP: Temperature.
     JS: List of 6 coupling constants.
     COLD-STEPS: Number of initial steps to discard (thermalization).
     RECORD-INTERVAL: How often to measure and save data (e.g., 10 = every 10th step).
     PRINT-INTERVAL: How often to print status to stdout.
     PRINT-SUMMARY: Boolean, if T print table to stdout.
   
   Returns (4 Lists):
     1. Energy History
     2. Magnetization Magnitude History
     3. Theta History
     4. Phi History"
  
  (initialize-spins *3d-spin-lattice*)
  
  (when (> cold-steps 0)
    (when print-summary (format t "Thermalizing for ~a steps...~%" cold-steps))
    (dotimes (i cold-steps)
      (metropolis-sweep *3d-spin-lattice* temp Js)))

  (let ((e-history nil)
        (mag-history nil)
        (theta-history nil)
        (phi-history nil))
    
    (when print-summary
      (format t "Starting Measurement Phase (T=~a, Steps=~a)...~%" temp steps)
      (format t "Step  | Energy   | Accept% | Mag |M| | Theta | Phi~%")
      (format t "----------------------------------------------------~%"))
    
    (dotimes (step steps)
      (let ((flips (metropolis-sweep *3d-spin-lattice* temp Js)))
        
        (let ((do-record (= (mod step record-interval) 0))
              (do-print  (and print-summary (= (mod step print-interval) 0))))
          
          (when (or do-record do-print)
            (let ((E (calculate-hamiltonian *3d-spin-lattice* Js)))
              (multiple-value-bind (avg-spin mag theta phi) 
                  (calculate-magnetization *3d-spin-lattice*)
                (declare (ignore avg-spin))
                
                (when do-record
                  (push E e-history)
                  (push mag mag-history)
                  (push theta theta-history)
                  (push phi phi-history))
                
                (when do-print
                  (let ((deg-theta (* theta (/ 180.0 pi)))
                        (deg-phi   (* phi   (/ 180.0 pi))))
                    (format t "~5d | ~8,2f | ~5,1f%  | ~5,3f | ~5,1f | ~5,1f~%" 
                            step E 
                            (* 100.0 (/ flips (array-total-size *3d-spin-lattice*)))
                            mag deg-theta deg-phi)))))))))
    
    (when print-summary (format t "Simulation Complete.~%"))
    
    (values (nreverse e-history)
            (nreverse mag-history)
            (nreverse theta-history)
            (nreverse phi-history))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;Ensemble Average & Physical Observables;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun analyze-series (data)
  "Calculates statistical properties of a dataset.
   
   Arguments:
     DATA: A list of numbers.
   
   Returns (2 Values):
     MEAN: Arithmetic average.
     VARIANCE: Population variance."
  (let ((sum 0.0) (sum-sq 0.0) (n (length data)))
    (if (= n 0) (values 0.0 0.0)
        (progn
          (dolist (x data)
            (incf sum x)
            (incf sum-sq (* x x)))
          (let ((mean (/ sum n)))
            (values mean (- (/ sum-sq n) (* mean mean))))))))

(defun ensemble-average (runs steps temp Js cold-steps record-interval)
  "Runs multiple independent simulations to calculate thermodynamic averages.
   
   Arguments:
     RUNS: Number of independent experiments.
     STEPS: Number of measurement steps per experiment.
     TEMP: Temperature.
     JS: Coupling constants.
     COLD-STEPS: Thermalization steps to discard.
     RECORD-INTERVAL: Sampling frequency.
   
   Returns (4 Values):
     1. <E>   : Ensemble Average Energy.
     2. Var(E): Variance of Energy (Specific Heat).
     3. <M>   : Ensemble Average Magnetization.
     4. Var(M): Variance of Magnetization (Susceptibility)."
  (let ((sum-E-mean 0.0) (sum-E-var 0.0)
        (sum-M-mean 0.0) (sum-M-var 0.0))
    
    (format t "Starting Ensemble Average (~a runs)...~%" runs)
    
    (dotimes (r runs)
      (multiple-value-bind (e-hist m-hist th ph) 
          (run-simulation steps temp Js cold-steps record-interval 1000 nil)
        (declare (ignore th ph))
        
        (multiple-value-bind (mean-e var-e) (analyze-series e-hist)
          (multiple-value-bind (mean-m var-m) (analyze-series m-hist)
            (incf sum-E-mean mean-e)
            (incf sum-E-var  var-e)
            (incf sum-M-mean mean-m)
            (incf sum-M-var  var-m)
            (format t "  Run ~2d/~d | <E>=~6,3f | <M>=~5,3f~%" (1+ r) runs mean-e mean-m)))))
    
    (format t "Ensemble Calculation Complete.~%")
    (values (/ sum-E-mean runs)
            (/ sum-E-var  runs)
            (/ sum-M-mean runs)
            (/ sum-M-var  runs))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;File I/O;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun multiple-temperature-run-data (filename start-T end-T step-T Js runs steps cold-steps record-interval)
  "Scans a range of temperatures and saves thermodynamic data to CSV.
   
   Arguments:
     FILENAME: Output file path (string).
     START-T, END-T, STEP-T: Temperature range parameters.
     JS: Coupling constants.
     RUNS: Ensemble size per temperature point.
     STEPS: Measurement steps per run.
     COLD-STEPS: Thermalization steps per run.
     RECORD-INTERVAL: Sampling frequency for measurements.
   
   Output:
     Writes to FILENAME. Columns: Temp, Energy, Cv, Magnetization, Susceptibility."
  
  (with-open-file (stream filename
                          :direction :output
                          :if-exists :supersede
                          :if-does-not-exist :create)
    
    (format stream "Temp,Energy,Cv,Magnetization,Susceptibility~%")
    (format t "Starting Temperature Scan (T=~,2f to ~,2f)...~%" start-T end-T)
    (do ((temp start-T (+ temp step-T)))
        ((> temp (+ end-T 0.0001)))
      (format t "~%--- Measuring T=~,3f ---~%" temp)
      
      (multiple-value-bind (avg-E var-E avg-M var-M)
          (ensemble-average runs steps temp Js cold-steps record-interval)
        
        (format stream "~,4f,~,5f,~,5f,~,5f,~,5f~%" 
                temp avg-E var-E avg-M var-M)
        (force-output stream))))
  (format t "~%Scan Complete. Data saved to ~a~%" filename))

(defun test-final-simulation ()
  "Runs a quick diagnostic test of the engine.
   Scans 3 points (T=0.5, 1.5, 2.5) to verify phase transition.
   Saves results to 'heisenberg_test.csv'.
   
   Output:
     Prints diagnostic to stdout and saves test CSV."
  
  (format t "==========================================~%")
  (format t "    TESTING HEISENBERG SIMULATION ENGINE  ~%")
  (format t "==========================================~%")
  
  (let ((filename "heisenberg_test.csv")
        (Js '(1.0 1.0 1.0 1.0 1.0 1.0)) ;; Standard Isotropic Interaction
        (start-T 0.5)
        (end-T   2.5)
        (step-T  0.25)
        (runs    2)    
        (steps   5000) 
        (cold    5000)
        (rec-int 50)) ;; User defined record interval
    
    (multiple-temperature-run-data 
     filename start-T end-T step-T Js runs steps cold rec-int)
    
    (format t "==========================================~%")
    (format t "Test Complete!~%")
    (format t "Output saved to: ~a~%" filename)
    (format t "Verify that 'Magnetization' drops from ~%~d (at T=0.5) to ~d (at T=2.5).~%" 
            "~1.0" "~0.0")
    (format t "==========================================~%")))

;;; HOW TO RUN:
;; (test-final-simulation)
