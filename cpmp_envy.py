import MainMELP as d
import time,sys,random,math
spp=0 #set partitioning
psize=1 #population size
maxloops=20 #max_loops_solution_not_improved
t=2000 #timelimit for searching
nf=0
n=0 # number of facilities, number of service areas
mthd=-1
distance_type=0
avg_dis=1.0
alpha=0
if len(sys.argv) >= 2:
    fn=sys.argv[1]
if len(sys.argv) >= 3:
    fn2=sys.argv[2]
if len(sys.argv) >= 4:
    nf=int(sys.argv[3])
if len(sys.argv) >= 5:
    avg_dis=float(sys.argv[4])
if len(sys.argv) >= 6:
    alpha=float(sys.argv[5])
if len(sys.argv) >= 7:
    psize=int(sys.argv[6])
if len(sys.argv) >= 8:
    maxloops=int(sys.argv[7])

d.location_problem=1 #1 for cpmp; 3 for pmp
d.distance_type=0
d.envy_service_objective=1 #0 for CPMP, 1 for CMELP or CMDELP
d.envy_objective_weight_of_travelcost=01 # 0 for CMELP; 1 for CMDELP 
d.envy_service_distance=avg_dis
d.envy_objective_weight=alpha
d.fixed_cost_obj=0
d.is_spp_modeling=0
d.adaptive_number_of_facilities=0
d.pop_dis_coeff=100000#10000#1000000 #100000 for FLP 1-5 for PDP
d.max_loops_solution_not_improved=maxloops #for search termination
d.ruin_oprators=[9]
d.initial_solution_method=1
d.distance_type=distance_type
if fn.find("ORLIB")>=0:
    d.read_bm_instance(fn) #read instances from or-lib(beasily),holmberge,yang...
elif fn.find(".txt")>=0:
    d.readfile(fn,fn2) #read self-defined instances
elif fn.find(".plc")>=0:
    d.read_bm_instance2(fn) #read tbed instances
#if fn2!="na": d.readdistance(fn2)

t0=time.time()
d.solution_similarity_limit=5.0 #max(10.0,100.0/n)
d.solver_message=0
d.max_num_facility=nf
d.mip_solver="gurobi"#"gurobi"#"cplex"
## MIP solver
#d.solver_message=01
#d.cpmp_mip_model(nf,100,1,t,0.0000001)

## matheuristic
d.cflpr_matheuristic(1000,100,0,psize,maxloops)

print "=========================Final results========================="
print "objective:",d.biobjective
print "facility cost",d.objective_fcost
print "transportation cost:",d.objective
print "srrvice overload", d.objective_overload
print "pool size",len(d.region_pool)
print "total time",time.time()-t0
print "facilities selected",[x for x in d.centersID if x>=0]
print "demand assignment:", d.node_groups
print "service area stat:"
for i in range(d.num_districts):
    if d.district_info[i][0]==0: continue
    print d.facilityCandidate[i],d.district_info[i], d.district_info[i][2]/d.district_info[i][1]
    

sol=d.sol_stat()
print "=========================Objectives========================="
for x in sol:
    print x,sol[x]

equality_measures=d.print_equality_measures()
print "-----------equality measures---------------"
for x in equality_measures:
    print  x,equality_measures[x]
print "CPU_time",time.time()-t0
