; particles
; particle dynamics by Runge-Kutta interpolation
;
; sample invocation
;
; by Otto Linnemann
; (C) 2011, GNU General Public Licence

(ns main
  (:refer-clojure)
  (:import (java.io BufferedWriter FileWriter))
  (:use [clojure.contrib.duck-streams :only (reader read-lines write-lines)])
  (:use [particles :only (solve gen_particle_states createLogger istates)] :reload)
  (:require clojure.contrib.string)
  )

(defn parse_initial_state "parse particles.ini file and generate initial state vector for solver"
  [filename]
  (let [ascii (slurp filename)
        lines (clojure.contrib.string/split-lines ascii)
        istates
        (map (fn [line]
               (let [line_elems (re-seq #"[-+[0-9\.]eE]+" line)
                     line_values (map #(Double/parseDouble %) line_elems)
                     pos  (take 3 line_values)
                     vel  (take 3 (drop 3 line_values))
                     mass (nth line_values 6)]
                 [pos vel mass]
                 ))
             (rest lines) ; skip first line (comment line)
             )]
    (let [vpos (map #(nth % 0) istates)
          vvel (map #(nth % 1) istates)
          vmass (map #(nth % 2) istates)]
      (gen_particle_states vpos vvel vmass))))


(defn write_logger_object [obj filename]
  (with-open [wtr (BufferedWriter. (FileWriter.	filename))]
    (.write wtr "plot \"-\", \"-\", \"-\"\n")
    (let [traj @(:trajectories_ref obj)]
      (dorun (map
               (fn [tj]  
                 (dorun 
                   (map (fn [pos_t] (.write wtr (apply format "%e %e %e %e\n" (map #(double %) pos_t))))
                        tj)) (.write wtr "e\n") 
                 ) 
               traj)))))


;(solve istates 10000 1e-5 @(def loggerObject (createLogger istates 1e-4)))
;(write_logger_object loggerObject "test.trj")


(def solar_init_states (parse_initial_state "data/solar.ini"))
(time (solve solar_init_states 3.1536e+08 50000 1e10 @(def loggerObject (createLogger solar_init_states 1e10))))
(write_logger_object loggerObject "solar.trj")
