(ns particles
  (:refer-clojure)
  (:use [clojure.contrib.generic 
         [math-functions :only [abs sqrt]]]
        )
  (:require clojure.stacktrace
            [clojure.contrib.generic.arithmetic :as ga])
  )

(def pos [[1 0 0] [0 -1 0] [-1 -1 0] [-1 0 0] [0 1 0]])
;(def vel (repeat (count pos) [0 0 0]))
;(def acc (repeat (count pos) [0 0 0]))

(defn pairs [pos] 
  (let [n (count pos)]
    (let [idx []]
      (for [s (range n) e (range (inc s) n)] 
        [s e]))))


; vector difference v1 - v2
(defn vdiff [v1 v2]
  (map #(- %1 %2) v1 v2))

; squared absolute value of vector sum ( v1^2, v2^2, ... )
(defn vabssquare [v]
  (reduce + (map #(* % %) v)))

; multiplies scalar s with vector v
(defn vsmult [v s]
  (map #(* s %) v))


(defn distance_func [v1 v2 G]
  (let [delta (vdiff v2 v1)
        squared_distance (vabssquare delta)
        distance (sqrt squared_distance)
        cubic_distance (* squared_distance distance)]
    (vsmult delta (/ G cubic_distance))
    ))

(defn gravity_pair_forces [pos G]
  (map #(distance_func (pos (% 0)) (pos (% 1)) G)
       (pairs pos)))

(defn gravity_acc [pos G]
  (let [pairs (pairs pos)
        pair_forces (map #(distance_func (pos (% 0)) (pos (% 1)) G) pairs)
        acc (vec (repeat (count pos) 0))]
; geht so nicht aber
    (do (map (fn [pair index]
           (let [source (pair 0) dest (pair 1)] 
             (assoc acc index "Hallo")
             ))  
         pairs (range (count pos))
         )
      acc)))

;aber so
;(vdiff
;  (filter #(== (% 0) pos_index) pairs) ;-> das abbilden auf pair forces
;  (filter #(== (% 1) pos_index) pairs)
;  )


; demo
(gravity_pair_forces pos 6.67429e-11) 

