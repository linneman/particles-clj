; particles
; particle dynamics by Runge-Kutta integration
;
; by Otto Linnemann
; (C) 2011, GNU General Public Licence

(ns particles
  (:refer-clojure :exclude [- + *])
  (:use [clojure.contrib.generic 
         [math-functions :only [sqrt]]
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

; squared absolute value
(defn sq_abs
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


;
; --- particle interaction ---
;

; gravity constant
(def G 6.67429e-11)

; computes all edges between indexes of position vector pos
; e.g. ([0 1] [0 2] [0 3] [0 4] [1 2] [1 3] [1 4] [2 3] [2 4] [3 4])
(defn pairs [pos] 
  (let [n (count pos)]
    (let [idx []]
      (for [s (range n) e (range (inc s) n)] 
        [s e]))))


; computes gravity function between two material points (mass not included)
(defn gravity_pair_force [v1 v2]
  (let [delta (- v2 v1)
        squared_distance (sq_abs delta)
        distance (sqrt squared_distance)
        cubic_distance (* squared_distance distance)]
    (* (/ G cubic_distance) delta )
    ))


; computes gravity function between all edges (material point pairs)
(defn gravity_pair_forces [pos]
  (map #(gravity_pair_force (pos (% 0)) (pos (% 1)))
       (pairs pos)))


; sums up a vector of vectors
(defn vsum [vi]
  (loop [res [0 0 0] v (vec vi)]
    (if (empty? v)
      res
      (recur (+ (vec res) (vec (first v))) (rest v)))))


; computes accellaration for all given positions and masses
(defn gravity_acc [pos masses]
  (let [pos (vec pos)
        masses (vec masses)
        pairs (pairs pos)
        pair_forces (map #(gravity_pair_force (pos (% 0)) (pos (% 1))) pairs)
        inter_force_pair (partition 2 (interleave pair_forces pairs))
        acc [0 0 0]]
    (map (fn [element]
           (let [
                 sources       (filter #(== ((second %) 0) element) inter_force_pair)
                 destinations  (filter #(== ((second %) 1) element) inter_force_pair)]
             (- 
               (vsum (map #(+ acc (* (masses ((second %) 1)) (first %))) sources))
               (vsum (map #(+ acc (* (masses ((second %) 0)) (first %))) destinations))))) 
         (range (count pos)))))
    



; --- Differential Equation Vectors for Solver  ---
;
; y: vector of the dynamic particle parameter input
; y[0] = pos.x
; y[1] = pos.y
; y[2] = pos.z
;
; y[3] = vel.x
; y[4] = vel.y
; y[5] = vel.z
;
; returns f: vector of the dynamic particle paramter output
; new velocity:
; f[0] = dy[0]/dt = e.g. y[3]
; f[1] = dy[1]/dt = e.g. y[4]
; f[2] = dy[2]/dt = e.g. y[5]
; new acceleration
; f[3] = dy[3]/dt = e.g. gravity_acc[0]
; f[4] = dy[4]/dt = e.g. gravity_acc[1]
; f[5] = dy[5]/dt = e.g. gravity_acc[2]
;


; generate vectors of y's
(defn gen_y_vectors [pos vel] 
  (partition 6 (flatten (interleave pos vel))))

; calculate vectors of f's
(defn dynamics_func_vectors [#^double t y_vec & {:keys [masses]}]
  (let [positions (map #(vec (take 3 %)) y_vec)
        vels (map #(drop 3 %) y_vec)
        accs (gravity_acc positions masses)]
   (partition 6 (flatten (interleave vels accs)))
  ))


; --- Runge Kutta ---
;


; --- demo ---
;
;initial conditions
(def pos [[1 0 0] [0 -1 0] [-1 -1 0] [-2 0 0] [0 1 4]])
(def vel (repeat (count pos) [0 0 0]))
(def masses (vec (repeat (count pos) 1.0)))

;(gravity_pair_forces pos) 
;(gravity_acc pos masses) 

(def initial_y_vectors (gen_y_vectors pos vel))
(def accs (dynamics_func_vectors 1.0 initial_y_vectors :masses masses))

initial_y_vectors
accs
