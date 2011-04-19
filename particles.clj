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


; computes all edges between indexes of position vector pos
; e.g. ([0 1] [0 2] [0 3] [0 4] [1 2] [1 3] [1 4] [2 3] [2 4] [3 4])
(defn pairs [pos] 
  (let [n (count pos)]
    (let [idx []]
      (for [s (range n) e (range (inc s) n)] 
        [s e]))))


; computes gravity function between two material points (mass not included)
(defn gravity_pair_force [v1 v2 G]
  (let [delta (- v2 v1)
        squared_distance (sq_abs delta)
        distance (sqrt squared_distance)
        cubic_distance (* squared_distance distance)]
    (* (/ G cubic_distance) delta )
    ))


; computes gravity function between all edges (material point pairs)
(defn gravity_pair_forces [pos G]
  (map #(gravity_pair_force (pos (% 0)) (pos (% 1)) G)
       (pairs pos)))


; sums up a vector of vectors
(defn vsum [vi]
  (loop [res [0 0 0] v (vec vi)]
    (if (empty? v)
      res
      (recur (+ (vec res) (vec (first v))) (rest v)))))


; computes accellaration for all given positions and masses
(defn gravity_acc [pos masses G]
  (let [pairs (pairs pos)
        pair_forces (map #(gravity_pair_force (pos (% 0)) (pos (% 1)) G) pairs)
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
    


;
; --- demo ---
;
(def G 6.67429e-11)
(def pos [[1 0 0] [0 -1 0] [-1 -1 0] [-2 0 0] [0 1 4]])
(def masses (vec (repeat (count pos) 1.0)))

;(def vel (repeat (count pos) [0 0 0]))
;(def acc (repeat (count pos) [0 0 0]))

;(gravity_pair_forces pos G) 
(gravity_acc pos masses G) 
