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


(def pos [[1 0 0] [0 -1 0] [-1 -1 0] [-2 0 0] [0 1 4]])
;(def vel (repeat (count pos) [0 0 0]))
;(def acc (repeat (count pos) [0 0 0]))



; vector algebra

(defmethod  -
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [a b]
  (map #(clojure.core/- %1 %2) a b))

(defmethod  -
  [clojure.lang.LazySeq clojure.lang.LazySeq]
  [a b]
  (map #(clojure.core/- %1 %2) a b))

(defmethod  +
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [a b]
  (map #(clojure.core/+ %1 %2) a b))

(defmethod  +
  [clojure.lang.LazySeq clojure.lang.LazySeq]
  [a b]
  (map #(clojure.core/+ %1 %2) a b))

(defmethod  *
  [java.lang.Number clojure.lang.PersistentVector]
  [s v]
  (map #(clojure.core/* s %) v))

(defmethod  *
  [java.lang.Number clojure.lang.LazySeq]
  [s v]
  (map #(clojure.core/* s %) v))

(defmethod  *
  [clojure.lang.PersistentVector clojure.lang.PersistentVector]
  [v1 v2]
  (reduce + (map #(* %1 %2) v1 v2)))


(defn sq_abs
  [x]
  (reduce + (map #(* % %) x)))


(defmulti abs class)

(defmethod abs
  java.lang.Number
  [x]
  (Math/abs x))

(defmethod abs
  clojure.lang.PersistentVector
  [x]
  (sq_abs x))



(defn pairs [pos] 
  (let [n (count pos)]
    (let [idx []]
      (for [s (range n) e (range (inc s) n)] 
        [s e]))))

(defn distance_func [v1 v2 G]
  (let [delta (- v2 v1)
        squared_distance (sq_abs delta)
        distance (sqrt squared_distance)
        cubic_distance (* squared_distance distance)]
    (* (/ G cubic_distance) delta )
    ))

(defn gravity_pair_forces [pos G]
  (map #(distance_func (pos (% 0)) (pos (% 1)) G)
       (pairs pos)))

; sums up a vector of vectors
(defn vsum [vi]
  (loop [res [0 0 0] v vi]
    (if (empty? v)
      res
      (recur (+ (vec res) (vec (first v))) (rest v)))))


(defn gravity_acc [pos G]
  (let [pairs (pairs pos)
        pair_forces (map #(distance_func (pos (% 0)) (pos (% 1)) G) pairs)
        inter_force_pair (partition 2 (interleave pair_forces pairs))
        acc (vec (repeat (count pos) 0))
        sum [0 0 0]]
    (map (fn [element]
           (let
             [source (filter #(== ((second %) 0) 2) inter_force_pair)
              dest (filter #(== ((second %) 1) 2) inter_force_pair)]
             (do (println element)
             (- 
               (vsum (map #(+ sum (vec (first %))) source))
               (vsum (map #(+ sum (vec (first %))) dest)))))) 
           (range (count pos)))))
    
; demo
;(gravity_pair_forces pos 6.67429e-11) 
(gravity_acc pos 10.0) 
