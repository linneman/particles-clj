; particles
; particle dynamics by Runge-Kutta integration
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

; particle data structure which provides state vector and
; constant characteristics e.g. its mass
(defstruct particles :state :mass)

; generate vectors of y state vector, see deqn
(defn gen_y_vectors [pos vel] 
  (partition 6 (flatten (interleave pos vel))))

; generates vector of particle states
(defn gen_particle_states [pos vel masses]
  (let [y_vectors  (gen_y_vectors pos vel)]
    (map #(struct particles (nth % 0) (nth % 1)) 
         (partition 2 (interleave y_vectors masses)))))


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
(defn gravity_deqn [t y p_reactio] 
  (let [acc (vsum
              (map #(* (:mass %)
                       (vec (gravity_pair_force
                              (vec (take 3 y))          ; actio position
                              (vec (take 3 (:state %))) ; reactio position
                              )))
                   p_reactio))]
    (vec (concat (vec (drop 3 y)) acc))
    ))


; --- Runge Kutta ---
;
(defn rkf45 [f t y h eps params]
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
        s (* 0.840896 (pow (/ (* eps h) (abs (- zj yj))) 0.25))
        ]
    {:yn zj :hn s} 
    )
  )


; --- demo ---
;
;initial conditions
(def pos [[1 0 -2] [0 -1 -1] [-1 -1 1] [-2 0 0] [0 1 4]])
(def vel (repeat (count pos) [0 0 0]))
(def masses (vec (repeat (count pos) 1.0)))
(def istates (gen_particle_states pos vel masses))

; test
(def actio (nth istates 0))
(def reactio 
  (keep-indexed (fn [idx item] (if (not (= idx 0)) item)) istates)
  )

; (defn rkf45 [f t y h eps params]

(defn mytest []
  ;(gravity_deqn  0.0 (:state actio) reactio)
  (rkf45 gravity_deqn 0.0 (:state actio) 300 1e-5 reactio)
  )




(defn indexed [col] (map vector (iterate inc 0) col))

(defn update [states t h eps]
  (pmap (fn [k]
         (let [actio (nth states k)
               reactio (keep-indexed (fn [idx item] (if (not (= idx k)) item)) states)
              ]
           (rkf45 gravity_deqn t (:state actio) h eps reactio)
           ))
       (range (count states))))


;(def initial_y_vectors (gen_y_vectors pos vel))
;initial_y_vectors

; examples for continuation
(def t 0)
(def h 300)
(def eps 1e-5)
(def yn_hn (update istates t h eps))
(def nstates (map (fn [a b] {:mass (:mass a) :state (:yn b)}) istates yn_hn))
(def trajectories (map #(vector (vec (take 3 (:state %)))) istates)) ; initialize with first position
(def lim_dist 1e-7)

(defn log [states]
 (do 
  (dorun (map #(println (take 3 (:state %))) states))
  (println "----------")))


(defn log2 [states] 
  (def trajectories
    (map (fn [a b] (conj b (take 3 (:state a)))) states trajectories)))

(defn log3 [states] 
  (def trajectories
    (map (fn [a b]
           (let [curr_pos (vec (take 3 (:state a)))
                 last_pos (vec (last b))
                 max_dist (apply max (map abs (- last_pos curr_pos)))]
             (if (> max_dist lim_dist) (conj b curr_pos) b)))
         states trajectories)))

(defn logt [states] 
  (def trajectories
    (map (fn [a b]
           (let [curr_pos (vec (take 3 (:state a)))
                 last_pos (vec (last b))
                 max_dist (apply max (map abs (- last_pos curr_pos)))]
                 b))
         states trajectories)))


(defn solve [tmax log_state_fn]
  (loop [t 0  states istates  h 300]
    (let [yn_hn (update nstates t h eps)
          updated_states (map (fn [a b] {:mass (:mass a) :state (:yn b)}) states yn_hn)
          h_rkf45 (apply max (map #(:hn %) yn_hn))
          [tn statesn] (if (> h h_rkf45)
                         (list t states)
                         (list (+ t h) (do (log_state_fn states) updated_states)))]
      (if (> t tmax)
        t
        (recur tn statesn (* 0.8 h_rkf45))
        ))))

;(solve 10000 log)
(solve 10000 log3)
