; particles
; particle dynamics by Runge-Kutta interpolation
;
; differential equation solver
;
; by Otto Linnemann
; (C) 2011, GNU General Public Licence

(ns particles
  (:refer-clojure :exclude [- + *])
  (:use [clojure.contrib.generic 
         [math-functions :only [sqrt pow]]
         [arithmetic :only [- + *]]
         ]
        )
  (:require clojure.stacktrace
            [clojure.contrib.generic.arithmetic :as ga])
  )

;
; --- vector algebra ---
;

; (v1 - v2)
(defmethod  -
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [a b]
  (vec (map #(clojure.core/- %1 %2) a b)))

; (v1 + v2)
(defmethod  +
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [a b]
  (vec (map #(clojure.core/+ %1 %2) a b)))

; (scalar + v)
(defmethod  +
  [java.lang.Number clojure.lang.PersistentVector]
  [s v]
  (vec (map #(clojure.core/+ s %) v)))

; (scalar * v)
(defmethod  *
  [java.lang.Number clojure.lang.PersistentVector]
  [s v]
  (vec (map #(clojure.core/* s %) v)))

; scalarproduct: (v1 * v2)
(defmethod  *
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [v1 v2]
  (reduce + (map #(* %1 %2) v1 v2)))

(defn sq_abs
  "squared absolute value of a vector"
  [x]
  (reduce + (map #(* % %) x)))

(defmulti abs class)

(defmethod abs
  java.lang.Number
  [x]
  (Math/abs x))

; absolute value of vector
(defmethod abs
  clojure.lang.PersistentVector
  [x]
  (sq_abs x))

(defn vsum
  "sums up a vector of vectors"
  [vi]
  (loop [res [0 0 0] v (vec vi)]
    (if (empty? v)
      res
      (recur (+ (vec res) (vec (first v))) (rest v)))))


;
; --- Particle Interaction and Differential Equation Vectors for Solver ---
;

; gravity constant
(def G 6.67429e-11)
(defn setG [val] (def G val))

; particle data structure which provides state vector and
; constant characteristics e.g. its mass
(defstruct particles :state :mass :t)


(defn gen_y_vectors
  "generate vectors of y state vector where y is the
   vector of the dynamic particle parameter input.
   refer also to function gravity_deqn.
   [0] = pos.x
   y[1] = pos.y
   y[2] = pos.z

   y[3] = vel.x
   y[4] = vel.y
   y[5] = vel.z"
  [pos vel] 
  (partition 6 (flatten (interleave pos vel))))


(defn gen_particle_states
  "generates a vector whose elements are pairs
  of particle state vectors y and corresponding
  masses. The input arguments are vectors of the
  same length for positions, velocities and
  masses."
  [pos vel masses]
  (let [y_vectors  (gen_y_vectors pos vel)]
    (map #(struct particles (nth % 0) (nth % 1) 0) 
         (partition 2 (interleave y_vectors masses)))))


(defn gravity_pair_force
  "computes gravity function between two material points (mass not included)"
  [v1 v2]
  (let [delta (- v2 v1)
        squared_distance (sq_abs delta)
        distance (sqrt squared_distance)
        cubic_distance (* squared_distance distance)]
    (* (/ G cubic_distance) delta )
    ))

(defn gravity_deqn 
  "Computes the differential dy/dt from a particle state y and
  the vector of all reacting particle states (positions). The
  functions returns the differential vector f (dynamic particle
  paramter output:

  new velocity:
  f[0] = dy[0]/dt = e.g. y[3]
  f[1] = dy[1]/dt = e.g. y[4]
  f[2] = dy[2]/dt = e.g. y[5]
  new acceleration
  f[3] = dy[3]/dt = e.g. gravity_acc[0]
  f[4] = dy[4]/dt = e.g. gravity_acc[1]
  f[5] = dy[5]/dt = e.g. gravity_acc[2]"
  [t y p_reactio] 
  (let [acc (vsum
              (map #(* (:mass %)
                       (vec (gravity_pair_force
                              (vec (take 3 y))          ; actio position
                              (vec (take 3 (:state %))) ; reactio position
                              )))
                   p_reactio))]
    (vec (concat (vec (drop 3 y)) acc))
    ))

;
; --- Runge Kutta ---
;
(defn rkf45
  "Runge Kutta forth/fifth order for the interpolation of 
  y(t+h) where h is the step size, t is the current time,
  y is the current state vector, f is the differential
  function, params are the arguments to that function and
  eps is the accepted error for the determination of the
  ideal step size hn. The ideal step size hn is calculated
  by taking into accout the given two interpolation orders.
  Both yn ( y(t+h) ) and hn are returned as structure
  elements."
  [f t y h eps params]
  (let [y (vec y)
        pf (fn [t y] (f t y params))
        k1 (* h (pf t y))
        k2 (* h (pf (+ t  (* 0.25 h)) (+ y (* 0.25 k1))))
        k3 (* h (pf (+ t  (* 0.375 h)) (+ y (* 0.09375 k1) (* 0.28125 k2))))
        k4 (* h (pf (+ t  (* 0.9230769230769231 h))
                    (+ y (* 0.8793809740555303 k1)
                       (* -3.277196176604461 k2)  (* 3.3208921256258535 k3))))
        k5 (* h (pf (+ t h)
                    (+ y (* 2.0324074074074074 k1) (* -8 k2)
                       (* 7.173489278752436 k3) (* -0.20589668615984405 k4))))
        k6 (* h (pf (+ t (* 0.5 h))
                    (+ y (* -0.2962962962962963 k1) (* 2 k2)
                       (* -1.3816764132553607 k3) (* 0.4529727095516569 k4)
                       (* 0.275 k5))))
        yj (+ y (* 0.11574074074074074 k1) (* 0.5489278752436647 k3)
              (* 0.5353313840155945 k4) (* -0.25 k5))
        zj (+ y (* 0.11851851851851852 k1) (* 0.5189863547758284 k3)
              (* 0.5061314903420167 k4) (* -0.18 k5) (* 0.03636363636363636 k6))
        tol (* eps (abs y))
        s (* 0.840896 (pow (/ (* tol h) (abs (- zj yj))) 0.25))
        ]
    {:yn zj :hn s} 
    )
  )


;
; --- Trajectory logging ---
;
(defprotocol Log
  "Protocol for Logging particle trajectories"
  (log [this states] "append current state information to log")) 

(defrecord Logger [trajectories_ref lim_dist]
  Log
  (log [this states] 
    (dosync (ref-set (:trajectories_ref this)
      (doall (map (fn [a b]
             (let [curr_pos (vec (take 3 (:state a)))
                   logged (vec (flatten [curr_pos (:t a)]))
                   last_pos (vec (last b))
                   max_dist (apply max (map abs (- last_pos curr_pos)))]
               (if (> max_dist (:lim_dist this)) (conj b logged) b)))
           states (deref (:trajectories_ref this)))))))
  )

(defn createLogger
  "creates a logger instance used for the solver. The function
  takes to arguments. istates correspond to the initial state
  vector and lim_dist provides the limit which must be crossed
  in one dimension for a new particles position for being logged." 
  [istates lim_dist]
  (let [trajectories (map #(vector (vec (flatten [(take 3 (:state %)) 0]))) istates)]
    (Logger. (ref trajectories) lim_dist)))



;
; --- Differential Equation Solver ---
;
(defn update
  "applies runge kutta function (rkf45) to all states in
  given state vector with step size h and error limit eps."
  [states t h eps]
  (pmap (fn [k]
         (let [actio (nth states k)
               reactio (keep-indexed (fn [idx item] (if (not (= idx k)) item)) states)
              ]
           (rkf45 gravity_deqn t (:state actio) h eps reactio)
           ))
       (range (count states))))

(defn solve
  "Solves initial value problem for istates. The step size is
  controlled automatically by means of Runge-Kutta interpolation
  in order to keep the local error smaller than eps and the
  interpolation process is stopped when the current time exceeds
  tmax. The transformation of the state vector is observed by an
  instance of Logger-object."
  ([istates tmax #^Log loggerObject]     (solve istates tmax 1e-6 loggerObject))
  ([istates tmax eps #^Log loggerObject] (solve istates tmax 1000000 eps loggerObject))
  ([istates tmax maxstep eps #^Log loggerObject]
  (loop [t 0  stepcnt 0  states istates  h 1.0]
    (let [yn_hn (update states t h eps)
          updated_states (map (fn [a b] {:mass (:mass a) :state (:yn b) :t t}) states yn_hn)
          h_rkf45 (apply max (map #(:hn %) yn_hn))
          log_state_fn (fn [states] (log loggerObject states))
          [tn stepcnt statesn] (if (> h h_rkf45)
                         (vector t stepcnt states)
                         (vector (+ t h) (inc stepcnt) (do (log_state_fn states) updated_states)))]
      (if (or (>= t tmax) (>= stepcnt maxstep))
        t
        (recur tn stepcnt (doall statesn) (* 0.8 h_rkf45))
        )))))


;
; --- Usage Illustration ---
;
;initial conditions
(def pos [[1 0 -2] [0 -1 -1] [-1 -1 1] [-2 0 0] [0 1 4]])
(def vel (repeat (count pos) [0 0 0]))
(def masses (vec (repeat (count pos) 1.0)))
(def istates (gen_particle_states pos vel masses))

; test
; (def actio (nth istates 0))
; (def reactio 
;   (keep-indexed (fn [idx item] (if (not (= idx 0)) item)) istates)
;   )

; (defn rkf45 [f t y h eps params]

; (defn mytest []
;   ;(gravity_deqn  0.0 (:state actio) reactio)
;   (rkf45 gravity_deqn 0.0 (:state actio) 300 1e-5 reactio)
;   )
;(solve 10000 log)
;(solve 10000 log3)
;(def myLogger (create_logger istates 1e-7))


; Invoke the solver
;(solve istates 10000 1e-5 @(def loggerObject (createLogger istates 1e-4)))

