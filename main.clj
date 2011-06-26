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

(defn parse_initial_state
  "parse particles.ini file and generate initial state vector for solver"
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


(defn parse_params
  "parse simulation parameter file"
  [filename]
  (let [ascii (slurp filename)
        extparfn (fn [parname]
                   (let [regexp (re-pattern
                     (format "(%s)(\\s+:=\\s+)([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)" parname))]
                     (Double/parseDouble (nth (re-find regexp ascii) 3))))
        kvpair (fn [parname] {(keyword parname) (extparfn parname)})]
    (particles/setG (extparfn "GRAVITY_CONST"))
    (merge
      (kvpair "GRAVITY_CONST")
      (kvpair "MAX_STEPS")
      (kvpair "MAX_TIME")
      (kvpair "EPS_REL")
      (kvpair "H_MAX")
      )
    ))


(defn write_logger_object [obj filename]
  "write simulation result from logger object obj to filename"
  (with-open [wtr (BufferedWriter. (FileWriter.	filename))]
    (.write wtr "plot \"-\", \"-\", \"-\"\n")
    (let [traj @(:trajectories_ref obj)]
      (dorun (map
               (fn [tj]  
                 (dorun 
                   (map (fn [pos_t] (.write wtr (apply format "%e %e %e %e\n"
                                                       (map #(double %) pos_t))))
                        tj)) (.write wtr "e\n") 
                 ) 
               traj)))))


;(solve istates 10000 1e-5 @(def loggerObject (createLogger istates 1e-4)))
;(write_logger_object loggerObject "test.trj")



;
; --- start the simulation ---
;

; load initial states
(def solar_init_states (parse_initial_state "data/solar.ini"))

; load simulation parameters
(def params (parse_params "data/solar.prm"))

; start solver
(time
  (solve
    solar_init_states
    (:MAX_TIME params)
    (:MAX_STEPS params)
    (:EPS_REL params)
    @(def loggerObject (createLogger solar_init_states
                                     (:H_MAX params))
       )))

; write results to output file
(write_logger_object loggerObject "solar.trj")

