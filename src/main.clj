; particles
; particle dynamics by Runge-Kutta interpolation
;
; main class for generation of stand-alone app
;
; by Otto Linnemann
; (C) 2011, GNU General Public Licence

(ns main
  (:gen-class)
  (:use [particles :only (solve createLogger)] :reload)
  (:use io))


; command line interface (leiningen)
(defn -main [& args]
  (if (not= (count args) 3)
    (println 
      "particle dynamics by Runge-Kutta interpolation
      invocation:  java -jar particles-0.1-standalone.jar  parameter-file  initials-states-file trajectory-output-file\n
      (C) 2011, GNU General Public Licence by Otto Linnemann")
    (let
      [pfilename (nth args 0) inifilename (nth args 1) trjfilename (nth args 2)]

      (dorun (println (format "parsing parameter file %s" pfilename)))
      (let [params (parse_params pfilename)]

        (dorun (println (format "parsing initial condition file %s" inifilename)))
        (let [initstates (parse_initial_state inifilename)]

          (dorun (println "start simulation ..."))

          ; start solver
          (time
            (solve
              initstates
              (:MAX_TIME params)
              (:MAX_STEPS params)
              (:EPS_REL params)
              @(def loggerObject (createLogger initstates
                                               (:H_MAX params))
                 )))

          ; write results to output file
          (write_logger_object loggerObject trjfilename)
          (System/exit 0)
          )))))


; usage illustration for testing
;(-main "data/solar.prm" "data/solar.ini" "data/solar.trj")
