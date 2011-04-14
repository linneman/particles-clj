(ns particles
  (:refer-clojure)
  (:use   (clojure.contrib.generic
            [math-functions :only [abs sqrt]]
            ))
  (:require clojure.stacktrace)
  )

(def pos [[1 0 0] [0 -1 0] [-1 -1 0] [-1 0 0] [0 1 0]])
;(def vel (repeat (count pos) [0 0 0]))
;(def acc (repeat (count pos) [0 0 0]))

(defn pairs [pos] 
  (let [n (count pos)]
    (let [idx []]
      (for [s (range n) e (range (inc s) n)] 
        [s e]))))

(defn deltas [pairs]
  (map #(for [coord (range 3)] 
          (- ((pos (% 0)) coord) ((pos (% 1)) coord))) 
       pairs))

(defn squared_distances [deltas]
  (map #(reduce + %)
       (map (fn dsquare [coord] (map #(* % %) coord))
            (vec deltas))))

(defn distances [squared_distances]
  (map #(sqrt %) (vec squared_distances)))


(defn cubic_distances [deltas]
  (let [sqd (squared_distances deltas) d (distances sqd)]
    (map #(* %1 %2) sqd d)))

(defn forces [pos gravity_const]
  (let [d (deltas (pairs pos)) cd (cubic_distances d)]
    (let [rec_cd (map #(/ 1.0 %) cd)]
      (map (fn [s v] (map #(* gravity_const s %) v)) rec_cd d)
      )))
 

; demo
(def d (deltas (pairs pos)))
(squared_distances d)
(distances (squared_distances d))
(cubic_distances d)

(forces pos 6.67429e-11) 

