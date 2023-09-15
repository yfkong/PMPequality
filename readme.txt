Computing enviroment:
Python 2.7
pypy 7.0 
pulp 2.4
Gurobi v9

PMP solver:
1 set d.envy_service_objective=0 in pmp_envy.py
2 run "pypy pmp_envy.py geo_zy.txt na 13"
command line: pypy pmp_envy.py datefilename na P

MELP solver:
1 set d.envy_service_objective=1, d.envy_objective_weight_of_travelcost=0,  in pmp_envy.py
2 run "pypy pmp_envy.py geo_zy.txt na 13 0.4 1.0"
command line: pypy pmp_envy.py datefilename na P d* weight_envy

MDELP solver:
1 set d.envy_service_objective=1, d.envy_objective_weight_of_travelcost=1, in pmp_envy.py
2 run "pypy pmp_envy.py geo_zy.txt na 13 0.4 1.0 5 100"
command line: pypy pmp_envy.py datefilename na P d* weight_envy num_solutions loops

For CPMP, CMELP and CMDELP, use "cpmp_envy.py" with the same parameters of pmp_envy.py.

 