# -*- coding: utf-8 -*-
## yfkong@henu.edu.cn, Apr.,2022
import sys,os,random,time,copy,math,tempfile
#mip solver
mip_solvers=[] #MIP solvers supported 
mip_solver=''  #MIP solver, "cplex", "cbc" or ""
mip_file_path=tempfile.gettempdir()
os.chdir(mip_file_path)  #used in arcgis
try:
    import pulp
    s=pulp.apis.GUROBI_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('gurobi')
    s=pulp.apis.cplex_api.CPLEX_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cplex')
    s=pulp.apis.COIN_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cbc')
except: 
    pass
if len(mip_solvers)>0: mip_solver=mip_solvers[0]

#constant
MAXNUMBER=1.0e+20
MINNUMBER=1.0e-10
#instance info
nodes=[]
nodes_std=[] #for homoDP only
weight_features=[] #for homoDP only
num_units=-1
nodedij=[]
nodedik=[]  #weighted cost from i to k, =nodedij*nodes[][3] 
nodendij=[] #network distance
node_neighbors=[]
facility_neighbors=[]
total_pop=0
avg_pop=0
total_supply=0
distance_type=0 #0 Euclidean, 1 Manhattan, 2 Geo 3 Network
all_units_as_candadate_locations=0
facilityCandidate=[] #attribute = 0,1,2...
facility_inclusion_list=[] #attribute = 1
facility_inclusion_list2=[] #attribute = 1

facilityCapacity=[]
facilityCost=[]
num_facilityCandidate=-1
num_districts=-1 # number of service areas/facilities
avg_dis_min=1.0
potential_facilities=[]
NearFacilityList=[]
nearCustomer=[]
nearCustomers=[]

#parameters for districting
location_problem=1
max_num_facility=999
adaptive_number_of_facilities=1
fixed_cost_obj=1
spatial_contiguity=0 # 0 no, 1 yes, 2 yes with multi-parts
pop_dis_coeff=10000.0 #used in the objective function
penalty_on_demand=10000.0

cflp_max_service_radius=10000.0
cflp_service_radius_preference=10000.0
cflp_min_service_coverage_percentage=100.0
envy_service_distance=10000.0
envy_objective_weight_of_travelcost=1 #pmp only
envy_objective_weight=0
envy_service_objective=0 
#0: no envy; 
#1: envy obj with weight; 
#2: envy constraint with CV; 
#3: envy with 0 weight on sum_dis, MELP
#4: envy(abs_dev) obj with weight; 
envy_coefficient_variation=0.5

#current solution
centersID=[]
node_groups=[]
district_info=[] #[[0,0,0.0] for x in range(num_districts)] # solution
objective_overload=0
obj_balance=MAXNUMBER
objective=MAXNUMBER
objective_fcost=0.0
biobjective=MAXNUMBER
objective_supply=0.0
objective_envy=0.0
objective_covered=0
objective_rmax_not_covered=0
objective_pentalty_on_overload=0.0
objective_pentalty_on_rmax=0.0
objective_pentalty_on_covering=0.0

given_solution=0 #reserved
all_solutions=[]

#best solution in each start
best_solution =[] # node_groups[:]
best_centersID=[]
best_biobjective=MAXNUMBER
best_objective=MAXNUMBER
best_objective_overload = MAXNUMBER
best_objective_fcost = 0.0
#global best solution 
#best_centers_global=[]
best_solution_global=[]
best_centersID_global=[]
best_biobjective_global = MAXNUMBER
best_objective_global = MAXNUMBER
best_objective_fcost_global = 0.0
best_overload_global = MAXNUMBER

#search statistics
time_check=0
time_check_edge_unit=0
time_spp=0.0
time_update_centers=0.0
time_op=[0.0 for x in range(10)]
time_ruin_recreate=[0.0 for x in range(10)]
time_location=[0.0 for x in range(10)]
time_pmp_re_location=0.0
time_Whitaker=0.0
time_repair=0
count_op=[0.0 for x in range(10)]
check_count=0
improved=0
move_count=0

#search histry
region_pool = []
pool_index=[]

#local search
assignment_operators_selected=[0,1] #0 one-unit move, 1 two-unit move, 2 three-unit move
location_operators_selected=[0,1,2,3,4] #0 swap, 1 drop, 2 add, 3 add+drop, 4 me
ruin_oprators=[0,1,2,3,4] #ruin0, ruin1, 9 mip assign
multi_start_count=6 #population size for GA, ILS, VNS, LIS+VND
initial_solution_method=0 #0 construction, 1 LP
assign_method=0 #not used
assign_or_Location_search_method=0
large_facility_cost=0
maxloops=1000
max_loops_solution_not_improved=-1
ruin_percentage=20
mainloop=0
solution_similarity_limit=10.0

#mip modelling for inititail solution, spp and full model
is_spp_modeling=1 #0 no spp modelling, 1 modelling at the end, 2 modelling in case of local optimum
linear_relaxation=0
spp_loops=400
solver_time_limit=7200 #used for mip modeling
solver_mipgap=0.000000000001
solver_message=0
heuristic_time_limit=300
seed =random.randint(0,10000)
random.seed(seed)

def arcpy_print(s):
    print s

#record a region in current solution
def update_region_pool(rid):
    global pool_index
    global time_spp
    global region_pool
    t=time.time()
    if is_spp_modeling<=0: return
    if centersID[rid]<0: return

    ulist=[x for x in  range(num_units) if node_groups[x]==rid]
    if ulist==[]:
        #print "empty area:",rid,node_groups
        return
    cost1=district_info[rid][2]
    cost2=sum(nodedik[x][rid] for x in ulist)
    if abs(cost1-cost2)>0.001: print rid,cost1,cost2
    obj=district_info[rid][2]+district_info[rid][4]*pop_dis_coeff
    idx=int(obj*100000)
    idx+=sum(x*(x+3) for x in ulist)
    if idx not in pool_index[rid]:
        pool_index[rid].append(idx)
        region_pool.append([ulist,district_info[rid][2],district_info[rid][1],district_info[rid][4],rid])
    time_spp+=time.time()-t
    return

#record all regions in current solution
def update_region_pool_all():
    if is_spp_modeling<=0:
        return
    for rid in range (num_districts):
        if centersID[rid]<0: continue
        update_region_pool(rid)

#update region information of the current solution
def update_district_info():
    global objective_overload
    global objective
    global biobjective
    global objective_fcost
    global district_info
    global move_count
    global obj_balance
    global centersID
    global objective_supply
    global avg_dis_min
    for k in range(num_districts):
        district_info[k][0] = 0
        district_info[k][1] = 0.0
        district_info[k][2] = 0.0
        district_info[k][3] = 0.0
        district_info[k][4] = 0.0
    for k in range(num_districts):
        if centersID[k]<0 and k in node_groups:
            arcpy_print("debug: a facility not selected but used: " + str(k))
            centersID[k]=facilityCandidate[k]
    for k in range(num_districts):
        if centersID[k]<0:
            continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        if len(ulist)==0:
            if location_problem==3: continue
            if adaptive_number_of_facilities==1 and k not in facility_inclusion_list:
                centersID[k]=-1
                continue
        district_info[k][0] = len(ulist)
        district_info[k][1] = sum(nodes[x][3] for x in ulist)
        district_info[k][2] = sum(nodedik[x][k] for x in ulist)
        if location_problem==9: district_info[k][2] = sum(nodedij[x][k] for x in ulist)
        district_info[k][3] = facilityCapacity[k] 
        district_info[k][4] = max(0.0,district_info[k][1]-facilityCapacity[k]) # -district_info[k][3]
        if location_problem==3: district_info[k][4]=0 #pmp
        #print centersID,node_groups
    bal=sum(x[4] for x in district_info)
    objective=sum([x[2] for x in district_info])
    objective_overload=bal
    #if objective/total_pop<avg_dis_min:
    #    avg_dis_min=objective/total_pop
    avg_dis_min=objective/total_pop
    biobjective=objective*envy_objective_weight_of_travelcost+objective_overload*avg_dis_min*pop_dis_coeff

    objective_supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x] >=0)
    #biobjective=objective+objective_overload*avg_dis_min*1000000
    #biobjective=bal2*avg_dis_min*1000000
    fcost=sum(facilityCost[x] for x in range(num_districts) if centersID[x] >=0)
    objective_fcost=fcost
    if fixed_cost_obj==1:
        biobjective+=fcost
    move_count+=1

def r_r_new_location(r): #add drop interchange swap/shift(add+drop) 
    global time_location
    global add_not_improved
    global location_tabu_list
    if len(location_operators_selected)==0: return -1,0
    nf=sum(1 for x in centersID if x>=0)
    location_operators=location_operators_selected[:]  #[0,1,2,3,4]#[0,1,3,4]
    #minregion=min(x[1] for x district_info if x[0]>0) 
    #if 0 in location_operators and (objective_supply-total_pop)/2 < minregion: location_operators.remove(0)
    idx=r
    if idx<0:
        idx=random.randint(0,len(location_operators)-1)
        idx=location_operators[idx]
    #if adaptive_number_of_facilities==0:
    #    if 1 in location_operators: location_operators.remove(1)
    #    if 2 in location_operators: location_operators.remove(2)
    sta=-1
    t=time.time()
    if idx==0:
        sta=r_r_location_swap_greedy()
        #time_location[0]+=time.time()-t
    if idx==1: 
        sta=r_r_location_drop_greedy() 
        #time_location[1]+=time.time()-t
    if idx==2: 
        sta=r_r_location_add_greedy()
        #time_location[2]+=time.time()-t
    if idx==3: 
        sta=r_r_location_add_greedy()
        RRT_local_search()
        sta=r_r_location_drop_greedy() 
        #time_location[3]+=time.time()-t
    if idx==4: 
        sta=r_r_reselect_location_univ()
        #time_location[4]+=time.time()-t
    if idx==5: 
        sta=r_r_reselect_location_pmp()
        #time_location[5]+=time.time()-t
    time_location[idx]+=time.time()-t
    return idx,sta

def r_r_location_drop_greedy(): #tested?
    global centersID
    global node_groups
    if adaptive_number_of_facilities==0: return 0
    if objective_supply<total_pop: 
        #print "drop objective_supply<total_pop: ", objective_supply 
        return 0
    dlist=[]
    clist=[x for x in range(num_districts) if centersID[x]>=0 and objective_supply-facilityCapacity[x]>=total_pop]
    if clist==[]: return 0
    for k in clist:
        ids=centersID[:]
        ids[k]=-1
        if locations_exist(ids)==1: continue
        if ids in location_tabu_list: continue
        ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x]==k]
        ulist.sort(key=lambda x:-x[1])
        savings=district_info[k][2]
        if fixed_cost_obj==1: savings+=facilityCost[k]
        caplist=[facilityCapacity[x]-district_info[x][1] for x in range(num_districts)]
        caplist[k]=0
        for x,y in ulist:
            assigned=0
            for r in NearFacilityList[x]:
                if centersID[r]<0: continue
                if r==k: continue
                if y<=caplist[r]:
                    savings-=nodedik[x][r]
                    caplist[r]-=y
                    assigned=1
                    break
            if assigned==0:
                for r in NearFacilityList[x]:
                    if centersID[r]<0: continue
                    if r==k: continue
                    #if caplist[r]<=0: continue
                    savings-=nodedik[x][r]
                    caplist[r]-=y
                    break
        surplus=sum(x for x in caplist if x <0) #negitive
        if surplus<-district_info[k][1]*0.05: continue
        savings+=surplus*avg_dis_min*5
        if savings<-biobjective*0.02: continue
        dlist.append([k,savings])
    if len(dlist)==0: return 0 #fail
    dlist.sort(key=lambda x:-x[1])
    if dlist[0][1]>0:
        while dlist[-1][1]<0: dlist.pop()
    r=random.random()
    idx=int(r*r*0.999*min(5,len(dlist)))
    k=dlist[idx][0]
    centersID[k]=-1
    
    ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x]==k]
    ulist.sort(key=lambda x:-x[1])
    caplist=[district_info[x][3]-district_info[x][1] for x in range(num_districts)]
    caplist[k]=0
    for x,y in ulist:
        assigned=0
        for r in NearFacilityList[x]:
            if centersID[r]<0: continue
            if y<=caplist[r]:
                caplist[r]-=y
                node_groups[x]=r
                assigned=1
                break
        if assigned==0:
            for r in NearFacilityList[x]:
                if centersID[r]<0: continue
                caplist[r]-=y
                node_groups[x]=r
                break
    #ulist=[x for x in range(num_units) if node_groups[x]==k]
    #for x in ulist:
    #    for k in NearFacilityList[x]:
    #        if centersID[k]>=0:
    #            node_groups[x]=k
    #            break
    update_district_info()
    #print "drop",dlist[idx][1],[biobjective,objective,objective_overload]
    return 1

def select_region(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    n=100*nf/num_units  #areas with 100 units
    if nf<=5: n=nf
    #if nf>=7: n=random.randint(nf/2+1,nf)
    if nf>=6 and nf<=11: n=random.randint(nf/2+1,9)
    if nf>=12 and nf<=16: 
        n=random.randint(7,10)
    if nf>=17: 
        n=random.randint(7,10)
    if n*num_units/nf<80: 
        n=min(10,80*nf/num_units)
    if location_problem==3: n=min(128,num_districts/10)
    u=seed
    if u<0: u=random.randint(0,num_units-1)
    r=node_groups[u]
    if objective_overload>0 and random.random()>0.2: #SAP
        rlist=[k for k in range(num_districts) if district_info[k][4]>0]
        random.shuffle(rlist)
        r=rlist[0]
        u=nearCustomer[r]
    clist=[r]
    for k in NearFacilityList[u]:
        if centersID[k]<0: continue
        if k==r: continue
        clist.append(k)
        if len(clist)==n: break
    return clist

def select_corridor_region():
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    n=100*nf/num_units  #areas with 100 units
    if nf<=5: 
        n=nf
        return [x for x in centersID if centersID[x]>=0]
    #if nf>=7: n=random.randint(nf/2+1,nf)
    if nf>=6 and nf<=11: n=random.randint(nf/2+1,9)
    if nf>=12 and nf<=16: 
        n=random.randint(7,10)
    if nf>=17: 
        n=random.randint(7,10)
    if n*num_units/nf<80: 
        n=min(10,80*nf/num_units)
    k_dis=[]
    
    for k in range(num_districts):
        if centersID[k]<0: continue
        if k not in node_groups:
            k_dis.append([k,MAXNUMBER])
            continue
        dis=max(nodedij[x][k] for x in range(num_units) if node_groups[x]==k)
        k_dis.append([k,dis])
    k_dis.sort(key=lambda x:-x[1])
    r=random.random()
    idx=int(r*r*len(k_dis)/2)
    k=k_dis[idx][0]
    clist=[k]
    while len(clist)<n:
        u=nearCustomer[clist[-1]]
        for k in NearFacilityList[u]:
            if centersID[k]<0: continue
            if k in clist: continue
            clist.append(k)
            break
    return clist

def select_region_cflpr(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    if nf<=10: n=random.randint(nf/2+1,nf)
    if nf>=11 and nf<=20:
         n=random.randint(8,10)
    if nf>=21: 
        n=random.randint(10,12)
    supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x]>=0)
    r=supply*1.0/total_pop
    if r>1.05: n=max(n, int(r/(r-1)))
    u=random.randint(0,num_units-1)
    r=node_groups[u]
    clist=[r]
    for k in NearFacilityList[u]:
        if centersID[k]<0: continue
        if k==r: continue
        clist.append(k)
        if len(clist)==n: break
    return clist

def select_region_adaptive(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    units=num_units/nf
    n=-1
    if units >= 50:
        if nf<=10: n=random.randint(nf/2+1,nf)
        if nf>=11 and nf<=20: n=random.randint(8,11)
        if nf>=21: n=random.randint(10,13)
    if units>=25 and units<=49:
        if nf<=10: n=random.randint(nf/2+1,nf)
        if nf>=11 and nf<=20: n=random.randint(9,12)
        if nf>=21: n=random.randint(11,14)
    if units <= 20:
        if nf<=10: n=random.randint(nf/2+1,nf)
        if nf>=11 and nf<=20: n=random.randint(10,13)
        if nf>=21: n=random.randint(12,15)
    n+= len(facility_inclusion_list)*n/nf
    if n>nf: n=nf
    u=seed
    if u<0 or u>=num_units: u=random.randint(0,num_units-1)
    r=node_groups[u]
    clist=[]
    for k in NearFacilityList[u]:
        if centersID[k]<0: continue
        clist.append(k)
        if len(clist)==n: break
    return clist

def get_cand_cflpr_locations(dlist,ulist,r):
    clist=[NearFacilityList[x][0] for x in ulist]
    clist=[x for x in clist if centersID[x]<0]
    clist=list(set(clist))
    random.shuffle(clist)
    n=int(len(dlist)*r)
    return clist[:n]

#update the best and the global best solution
#if the current solution is better than the best
def update_best_solution():
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_fcost
    global best_overload
    global best_objective_overload
    global best_centersID_global
    global best_solution_global
    global best_biobjective_global
    global best_objective_global
    global best_objective_fcost_global
    global best_overload_global    
    global improved_loop
    global improved
    global avg_dis_min
    improve =0
    if location_problem==1 and adaptive_number_of_facilities==0:
        nf=sum(1 for x in centersID if x>=0)
        if nf!=max_num_facility: return 0
    biobj=biobjective
    biobj_best=best_biobjective
    if biobj<=biobj_best:
        best_biobjective=biobj
        best_objective = objective
        best_objective_fcost=objective_fcost
        best_objective_overload = objective_overload
        best_solution = node_groups[:]
        best_centersID=centersID[:]
        improved_loop=mainloop
        improve =1
        improved+=1
    if biobj<best_biobjective_global:
        best_biobjective_global=biobj
        best_centersID_global=centersID[:]
        best_overload_global = objective_overload
        best_solution_global =node_groups[:]
        best_objective_global = objective
        best_objective_fcost_global=objective_fcost
        avg_dis_min=biobj/total_pop
    return improve

#return the neighor regions of unit nid
def neighbor_regions(nid):
    rid=node_groups[nid]
    knn=4
    rlist=[]
    for k in NearFacilityList[nid]:
            if k==rid: continue
            if centersID[k]<0: continue
            rlist.append(k)
            if len(rlist)>=knn: break 
    return rlist

def one_unit_move(): #FA
    global node_groups
    global time_op
    global count_op
    t=time.time()
    ulist=range(num_units)
    difflist=[]
    for x in ulist:
        rlist=[r for r in NearFacilityList[x] if centersID[r]>=0]
        k=node_groups[x]
        difflist.append([x,nodedik[x][rlist[0]]-nodedik[x][k]])
    difflist.sort(key=lambda x:x[1])
    ulist=[x[0] for x in difflist]
    find_best_solution=0
    for nid in ulist:
        rid=node_groups[nid]
        demand=nodes[nid][3]
        klist=neighbor_regions(nid)
        if len(klist)>1:
            klist2=[x for x in NearFacilityList[nid] if x in klist]
            klist=klist2
        for k in klist:
            #if k==rid: print "k==rid", nid,rid,neighbor_regions(nid)
            if demand+district_info[k][1]>facilityCapacity[k]: continue
            savings=nodedik[nid][rid]-nodedik[nid][k]
            overload_old=district_info[rid][4]+district_info[k][4]
            overload_new= max(0,district_info[rid][1]-nodes[nid][3]-facilityCapacity[rid]) +max(0,district_info[k][1]+nodes[nid][3]-facilityCapacity[k])
            savings+= (overload_old-overload_new)*avg_dis_min*pop_dis_coeff
            if envy_service_objective==1:
                if nodedij[nid][rid]>envy_service_distance:
                    savings+= nodes[nid][3]*envy_objective_weight*(nodedij[nid][rid]-envy_service_distance)*(nodedij[nid][rid]-envy_service_distance)
                if nodedij[nid][k]>envy_service_distance:
                    savings-= nodes[nid][3]*envy_objective_weight*(nodedij[nid][k]-envy_service_distance)*(nodedij[nid][k]-envy_service_distance)
            if envy_service_objective==4:
                if nodedij[nid][rid]>envy_service_distance:
                    savings+= nodes[nid][3]*envy_objective_weight*(nodedij[nid][rid]-envy_service_distance)
                if nodedij[nid][k]>envy_service_distance:
                    savings-= nodes[nid][3]*envy_objective_weight*(nodedij[nid][k]-envy_service_distance)
            if savings<=0: continue
            node_groups[nid]=k
            count_op[0]+=1
            #print "OUM",best_savings,biobjective,objective,objective_overload,
            update_district_info()
            #print "->",biobjective,objective,objective_overload
            find_best_solution += 1#update_best_solution()
            break
    time_op[0]+=time.time()-t
    return find_best_solution

def pmp_VND_search():
    global node_groups
    global time_op
    global count_op
    improved = 0
    t=time.time()
    while 1:
        improve = 0
        ulist=find_edge_units()
        for u in ulist:
            k=node_groups[u]
            for u2 in node_neighbors[u]:
                if node_groups[u2]==k: continue
                k2=node_groups[u2]
                if nodedik[u][k2]>=nodedik[u][k]: continue
                node_groups[u] = k2
                improve+=1
                improved+=1
        if improve == 0: break
    time_op[0]+=time.time()-t
    #print "improved",improved,
    update_district_info()
    return improved

def two_unit_move():
    global node_groups
    global time_op
    global count_op
    t=time.time()
    find_best_solution=0
    improve = 0
    ulist=range(num_units)
    difflist=[]
    for x in ulist:
        rlist=[r for r in NearFacilityList[x] if centersID[r]>=0]
        k=node_groups[x]
        difflist.append([x,nodedij[x][rlist[0]]-nodedij[x][k]])
    if len(difflist)==0:
        print "debug: two_unit_move(): len(difflist)==0"
    difflist.sort(key=lambda x:x[1])
    ulist=[x[0] for x in difflist]
    nu=len(ulist)/2
    fulist=ulist[:nu]
    movelist =[]
    for n1 in ulist:
        k1 = node_groups[n1]
        rlist1=neighbor_regions(n1)
        success_move=0
        ulist2=[x for x in ulist if node_groups[x] in rlist1]
        random.shuffle(ulist2)
        for n2 in ulist2:
            #if n1 == n2: continue
            k2=node_groups[n2]
            #if k2 not in rlist1: continue
            #if district_info[k2][1]+nodes[n1][3]-nodes[n2][3]-facilityCapacity[k2] >district_info[k2][4]: continue
            for k3 in neighbor_regions(n2):
                #if k3!=k1 and district_info[k3][1]+nodes[n2][3]>facilityCapacity[k3]: continue
                #if k3==k1 and district_info[k3][1]-nodes[n1][3]+nodes[n2][3]-facilityCapacity[k3]>district_info[k3][4]: continue
                savings=nodedik[n1][k1]+nodedik[n2][k2]
                savings-=nodedik[n1][k2]+nodedik[n2][k3]
                if envy_service_objective==1:
                    if nodedij[n1][k1]>envy_service_distance:
                        savings+= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k1]-envy_service_distance)*(nodedij[n1][k1]-envy_service_distance)
                    if nodedij[n2][k2]>envy_service_distance:
                        savings+= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k2]-envy_service_distance)*(nodedij[n2][k2]-envy_service_distance)
                    if nodedij[n1][k2]>envy_service_distance:
                        savings-= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k2]-envy_service_distance)*(nodedij[n1][k2]-envy_service_distance)
                    if nodedij[n2][k3]>envy_service_distance:
                        savings-= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k3]-envy_service_distance)*(nodedij[n2][k3]-envy_service_distance)
                if envy_service_objective==4:
                    if nodedij[n1][k1]>envy_service_distance:
                        savings+= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k1]-envy_service_distance)
                    if nodedij[n2][k2]>envy_service_distance:
                        savings+= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k2]-envy_service_distance)
                    if nodedij[n1][k2]>envy_service_distance:
                        savings-= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k2]-envy_service_distance)
                    if nodedij[n2][k3]>envy_service_distance:
                        savings-= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k3]-envy_service_distance)

                savings2=0
                if k3!=k1:
                    slist=[district_info[k1][1]-nodes[n1][3],district_info[k2][1]+nodes[n1][3]-nodes[n2][3],district_info[k3][1]+nodes[n2][3]]
                    savings2=district_info[k1][4]+district_info[k2][4]+district_info[k3][4]
                    savings2-=max(0,slist[0]-facilityCapacity[k1])+max(0,slist[1]-facilityCapacity[k2])+max(0,slist[2]-facilityCapacity[k3])
                if k3==k1:
                    slist=[district_info[k1][1]-nodes[n1][3]+nodes[n2][3],district_info[k2][1]+nodes[n1][3]-nodes[n2][3]]
                    savings2=district_info[k1][4]+district_info[k2][4]
                    savings2-=max(0,slist[0]-facilityCapacity[k1])+max(0,slist[1]-facilityCapacity[k2])
                #if district_info[k1][4]==0 and savings<0: continue
                if savings2<0: continue
                if savings2==0 and savings<0: continue
                #if savings<=0: continue
                count_op[1]+=1
                sp=objective_supply
                node_groups[n1] = k2
                node_groups[n2] = k3
                obj=biobjective
                #print "TUM",biobjective,objective,objective_overload,
                sp=objective_supply
                update_district_info()
                #print "->",biobjective,objective,objective_overload,":",savings,savings2,biobjective-obj
                find_best_solution += 1#update_best_solution()
                success_move=1
                break
            if success_move==1: break            
    time_op[1]+=time.time()-t
    return find_best_solution


# local search
def RRT_local_search():
    #if location_problem==3:
    #    pmp_VND_search()
    #    return
    global improved
    #global node_neighbors
    improved=0
    #operators=[0,1,2,3,4]
    operators=assignment_operators_selected[:]
    #if op_random == 1 and random.random()>0.5:
    #    random.shuffle(operators)
    for op in operators:
        if op == 0:
            one_unit_move()
            #update_region_pool_all()
        if op == 1:
            two_unit_move()
            #update_region_pool_all()
        if op == 2:
            #if random.random()<0.1: three_unit_move()
            three_unit_move()
            #update_region_pool_all()
    return

#local search with operator op
def vnd_op_search(op):
    #global node_neighbors
    #for x in node_neighbors:
    #    random.shuffle(x)
    if op == 0:
        one_unit_move()
    if op == 1:
        two_unit_move()
    #if 2 in assignment_operators_selected: print "VND",op, biobjective
    return

#VND search
def VND_local_search():
    #if location_problem==3:
    #    pmp_VND_search()
    #    return
    global improved
    improved=0
    #operators=[0,1,2,3,4]
    operators=assignment_operators_selected[:]    
    if len(operators)==0: return
    #if op_random == 1:
    #if op_random == 1 and random.random()>0.5:
    #    random.shuffle(operators)
    obj=biobjective
    while 1:
        vnd_op_search(operators[0])
        if biobjective < obj-0.00001:
            obj=biobjective
            continue
        if len(operators)==1:break
        vnd_op_search(operators[1])
        if biobjective < obj-0.00001:
            obj=biobjective
            continue
        if len(operators)==2:break

        vnd_op_search(operators[2])
        if biobjective < obj-0.00001:
            obj=biobjective
            continue
        if len(operators)==3:break
        vnd_op_search(operators[3])
        if biobjective < obj-0.00001:
            obj=biobjective
            continue
        if len(operators)==4:break
        vnd_op_search(operators[4])
        if biobjective < obj-0.00001:
            obj=biobjective
            continue
        break
    return

def read_bm_instance(f1):
    global num_units
    global total_pop
    global total_supply
    global nodes
    global node_neighbors
    global nodedij
    global nodedik
    global centersID
    global facilityCandidate
    global facilityCapacity
    global facilityCost
    global num_facilityCandidate
    global num_districts
    global district_info
    global avg_dis_min
    global potential_facilities
    node =[0,0,0,0,0,0,0,0,0,0]
    #school_nodes = []
    nodes = []
    f = open(f1)
    line = f.readline() #I,J
    line=line[0:-1]
    items = line.split(' ')
    idx=0
    for item in items:
        if item=="":
            continue
        if idx==0:
            num_districts=int(item)
            idx+=1
        else:
            num_units=int(item)
            idx+=1
    facilityCandidate=[x for x in range(num_districts)]
    facilityCapacity=[0.0 for x in range(num_districts)]
    facilityCost=[0.0 for x in range(num_districts)]
    nodes=[node[:] for x in range(num_units) ]
    nodedik=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
    nodedij=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
    arcpy_print("M,N: "+str(num_districts)+" "+str(num_units))
    for i in range(num_districts):
        line = f.readline()
        line=line[0:-1]
        items = line.split(' ')
        idx=0
        for item in items:
            if item=="":
                continue
            if idx==0:
                facilityCapacity[i]=float(item)
                idx+=1
            else:
                facilityCost[i]=float(item)
                idx+=1
    idx=0
    line = f.readline()
    while line: # for i in range((num_units+1)/10):
        line=line[0:-1]
        items = line.split(' ')
        for item in items:
            if item=="":
                continue
            if idx<num_units:
                nodes[idx][3]=float(item)
                idx+=1
            else:
                #i=(idx-num_units)/num_districts
                #k=(idx-num_units)%num_districts
                i=(idx-num_units)%num_units
                k=(idx-num_units)/num_units
                nodedik[i][k]=float(item)
                nodedij[i][k]=float(item)/nodes[i][3]
                idx+=1
        line = f.readline()    
    f.close()

    centersID=facilityCandidate[:]
    total_pop=sum(x[3] for x in nodes)
    total_supply=sum(facilityCapacity)
    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    avg_dis_min=1.0
    create_facility_neighbors()
    find_NearFacilityList(num_districts)
    #print NearFacilityList
    find_near_customer()
    create_node_neighbors()
    #find_nearFacilityFacility()
    potential_facilities=[x for x in range(num_districts)]
    s="total demand: "+ str(total_pop)
    arcpy_print(s)
    s="total supply: "+str(total_supply)
    arcpy_print(s)

def read_bm_instance2(f1):
    global num_units
    global total_pop
    global total_supply
    global nodes
    global node_neighbors
    global nodedij
    global nodedik
    global centersID
    global facilityCandidate
    global facilityCapacity
    global facilityCost
    global num_facilityCandidate
    global num_districts
    global district_info
    global avg_dis_min
    global potential_facilities
    node =[0,0,0,0,0,0,0,0,0,0]
    #school_nodes = []
    nodes = []
    f = open(f1)
    line = f.readline() #I,J
    line=line[0:-1]
    items = line.split(' ')
    if line.find("\t"): items = line.split('\t')
    idx=0
    for item in items:
        if item=="":
            continue
        if idx==0:
            num_units=int(item)
            idx+=1
        else:
            num_districts=int(item)
            idx+=1
    facilityCandidate=[x for x in range(num_districts)]
    facilityCapacity=[0.0 for x in range(num_districts)]
    facilityCost=[0.0 for x in range(num_districts)]
    nodes=[node[:] for x in range(num_units) ]
    nodedik=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
    nodedij=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
    arcpy_print("M,N: "+str(num_districts)+" "+str(num_units))

    idx=0
    while 1:
        line = f.readline()
        line=line[0:-1]
        if line=="": continue
        items = line.split(' ')
        if line.find("\t")>=0: items = line.split('\t')
        for item in items:
            if item=="":
                continue
            nodes[idx][3]=float(item)
            idx+=1
        if idx==num_units: break

    fidx=0
    while 1:
        line = f.readline()
        line=line[0:-1]
        if line=="": continue
        items = line.split(' ')
        if line.find("\t")>=0: items = line.split('\t')
        for item in items:
            if item=="":
                continue
            facilityCapacity[fidx]=float(item)
            fidx+=1
        if fidx==num_districts: break
    fidx=0
    while 1:
        line = f.readline()
        line=line[0:-1]
        if line=="": continue
        items = line.split(' ')
        if line.find("\t")>=0: items = line.split('\t')
        for item in items:
            if item=="":
                continue
            facilityCost[fidx]=float(item)
            fidx+=1
        if fidx==num_districts: break
    idx=0
    while 1:
        line = f.readline()
        line=line[0:-1]
        if line=="": continue
        items = line.split(' ')
        if line.find("\t")>=0: items = line.split('\t')
        for item in items:
            if item=="":
                continue
            k=idx/num_units
            i=idx%num_units
            idx+=1
            nodedik[i][k]=float(item)*nodes[i][3]
            nodedij[i][k]=float(item)
        if idx==num_units*num_districts: break
    f.close()
    centersID=facilityCandidate[:]
    total_pop=sum(x[3] for x in nodes)
    total_supply=sum(facilityCapacity)
    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    avg_dis_min=1.0
    create_facility_neighbors()
    find_NearFacilityList(num_districts)
    find_near_customer()
    #find_nearFacilityFacility()
    create_node_neighbors()
    potential_facilities=[x for x in range(num_districts)]
    s="total demand: "+ str(total_pop)
    arcpy_print(s)
    s="total supply: "+str(total_supply)
    arcpy_print(s)

def create_facility_neighbors():
    return 
    global facility_neighbors
    mindij=[[MAXNUMBER for x in range(num_districts)] for y in range(num_districts)]
    for i in range(num_districts):
        for j in range(num_districts):
            if j<=i: continue
            dlist=[nodedij[x][i]-nodedij[x][j] for x in range(num_units)]
            d=sum(x*x for x in dlist)
            mindij[i][j]=d
            mindij[j][i]=d
    facility_neighbors = [[]for x in range(num_districts)]
    for i in range(num_districts):
        dlist=[[x, mindij[i][x]] for x in range(num_districts) ]
        dlist.sort(key=lambda x:x[1])
        nghrs=[x[0] for x in dlist]
        facility_neighbors[i]=nghrs[:]

def create_node_neighbors():
    global node_neighbors
    #rlist=[x for x in range(num_districts)]
    mindij=[[MAXNUMBER for x in range(num_units)] for y in range(num_units)]
    for i in range(num_units):
        for j in range(num_units):
            if j<=i: continue
            dlist=[nodedij[i][x]-nodedij[j][x] for x in range(num_districts)]
            d=sum(x*x for x in dlist)
            mindij[i][j]=d
            mindij[j][i]=d
    node_neighbors = [[]for x in range(num_units)]
    for i in range(num_units):
        dlist=[[x, mindij[i][x]] for x in range(num_units) ]
        dlist.sort(key=lambda x:x[1])
        nn=8
        if nn>num_units: nn=num_units
        nghrs=[dlist[x][0] for x in range(nn)]
        random.shuffle(nghrs) #debug
        node_neighbors[i]=nghrs[:]

#read instance file, f1:unit info, f2: connectivity info
def readfile(f1,f2):
    global num_units
    global total_pop
    global total_supply
    global nodes
    global node_neighbors
    global nodedij
    global nodedik
    global centersID
    global facilityCandidate
    global facilityCapacity
    global facilityCost
    global num_facilityCandidate
    global num_districts
    global district_info
    global avg_dis_min
    global potential_facilities
    global facility_inclusion_list
    node =[0,0,0,0,0,0,0,0,0,0]
    #school_nodes = []
    nodes = []
    #nodes.append(node)
    print "reading nodes ...",
    f = open(f1)
    line = f.readline()  #OID    pop    PointX    PointY    fcadidature    fcost    fcapacity
    line = f.readline()
    nodeidx=0
    while line:
        line=line[0:-1]
        #print line
        items = line.split(',')
        if len(items)<=2:
            items = line.split('\t')
        unit=[nodeidx, float(items[2]), float(items[3]), int(items[1]),int(items[0]),int(items[6]),int(items[4]),float(items[5])]
        nodes.append(unit)
        nodeidx+=1
        #nodes.append([int(items[1]), float(items[8]), float(items[9]), int(items[5]), int(items[6]), int(items[7]),int(items[12]),int(items[13])])
        line = f.readline()
    f.close()
    num_units=len(nodes)
    total_pop=sum(x[3] for x in nodes)
    ##noprint num_units,"units"
    ##noprint "reading connectivity ...",
    ###id1,id2#####
    node_neighbors = [[]for x in range(len(nodes))]
    if f2!="na":
        #connectivity=[[0 for x in range(len(nodes))] for y in range(len(nodes))]
        f = open(f2)
        line = f.readline()
        line = f.readline()
        links=0
        while line:
            items = line.split(',')
            if len(items)<=2:
                items = line.split('\t')
            if int (items[1]) != int (items[2]):
                id1=int (items[1])
                id2=int (items[2])
                idx1=-1
                idx2=-1
                for i in range(num_units):
                    if nodes[i][4]==id1:
                        idx1=i
                    if nodes[i][4]==id2:
                        idx2=i
                    if idx1>=0 and idx2>0:
                        break
                node_neighbors[idx1].append(idx2)
                node_neighbors[idx2].append(idx1)
                #connectivity[idx1][idx2]=1
                #connectivity[idx2][idx1]=1
                links+=1
            line = f.readline()
        f.close()
    ##noprint links,"links"
    num_units=len(nodes)
    facilityCandidate=[]
    facilityCapacity=[]
    facilityCost=[]
    centersID=[]
    print "all data are read! "
    for i in range(num_units):
        if nodes[i][5]>0 or all_units_as_candadate_locations==1:
            facilityCandidate.append(i)
            facilityCapacity.append(nodes[i][5])
            facilityCost.append(nodes[i][7])
            centersID.append(i)
    num_facilityCandidate=len(facilityCandidate)
    num_districts=len(facilityCandidate)
    #facilityCandidate.sort()
    total_supply=sum(facilityCapacity)
    centersID=facilityCandidate[:]
    facility_inclusion_list=[]
    for k in range(num_districts):
        u=facilityCandidate[k]
        if nodes[u][6]==1:
            facility_inclusion_list.append(k)
    print "existing facilities", facility_inclusion_list
    nodedij=[[MAXNUMBER for x in range(num_districts)] for y in range(num_units)]
    max_dij=0.0
    print "craeating distance matrix......"
    for i in range(num_units):
        for k in range(num_districts):
            j=facilityCandidate[k]
            d=0.0
            if distance_type==0:
                d2=pow(nodes[i][1]-nodes[j][1],2)
                d2+=pow(nodes[i][2]-nodes[j][2],2)
                d=pow(d2,0.5)/1000
            if distance_type==1:
                d=abs(nodes[i][1]-nodes[j][1])
                d+=abs(nodes[i][2]-nodes[j][2])
                d/=1000
            if distance_type==2: #geo
                d=abs(nodes[i][1]-nodes[j][1])
                d+=abs(nodes[i][2]-nodes[j][2])
                d/=1000
            nodedij[i][k]=d
            if d>max_dij:
                max_dij=d
    ''' 
    if len(node_neighbors)>0:
        for i in range(len(nodes)):
            for j in range(len(nodes)):
                if j<=i:
                    continue
                if connectivity[i][j]==1:
                    node_neighbors[i].append(j)
                    node_neighbors[j].append(i)'''

    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    dis=0.0
    for i in range(num_units):
        d=min(nodedij[i])
        dis+=d*nodes[i][3]
    avg_dis_min=dis/total_pop
    #weighted cost from i to k
    
    nodedik=[[nodedij[y][x]*nodes[y][3] for x in range(num_districts)] for y in range(num_units)]
    find_NearFacilityList(num_districts)
    print "find_near_customer()..."
    find_near_customer()
    #find_nearFacilityFacility()
    print "create_facility_neighbors()..."
    create_facility_neighbors()
    potential_facilities=[x for x in range(num_districts)]
    s="M N: "+str(num_districts)+" "+str(num_units)
    arcpy_print(s)
    s="Total demand & supply: "+str(total_pop)+" "+str(total_supply)+" "+str(sum(facilityCost))
    arcpy_print(s)

def find_nearFacilityFacility():
    global nearFacilityFacility
    nearFacilityFacility=[[] for x in range(num_districts)]
    dkk=[sum(nodedik[x][k]*nodedik[x][k] for x in range(num_units)) for k in range(num_districts)]
    #dkk.sort(key=lambda x:x[1])
    #dlist=[x[0] for x in dkk]
    for k in range(num_districts):
        d=dkk[k]
        dk=[[x,dkk[x]-d] for x in range(num_districts)]
        dk.sort(key=lambda x:x[1])
        del dk[0]
        nearFacilityFacility[k]=[x[0] for x in dk]
   
def find_near_customer():
    global nearCustomer
    global nearCustomers
    nearCustomer=[-1 for x in range(num_districts)]
    nearCustomers=[[] for x in range(num_districts)]
    dis=[]
    for k in range(num_districts):
        dis=[ [x,nodedij[x][k]] for x in range(num_units)]
        dis.sort(key=lambda x: x[1])
        nearCustomer[k]=dis[0][0]
        nearCustomers[k]=[x[0] for x in dis]
       
def initialize_instance():
    global num_units
    global num_districts
    global num_facilityCandidate
    global centersID
    global node_groups
    global facilityCost
    global facilityCandidate
    global facilityCapacity
    global nodedik
    global avg_pop
    global total_pop
    global avg_dis_min
    global total_supply
    global fixed_cost_obj
    global max_num_facility
    global adaptive_number_of_facilities
    #solution obj 
    global district_info
    global objective_overload
    global objective
    global biobjective
    global all_solutions
    #best solution 
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_overload
    #global best solution 
    global best_solution_global
    global best_centersID_global
    global best_biobjective_global
    global best_objective_global
    global best_overload_global
    global potential_facilities
    global max_exclusion_list
    num_districts=len(facilityCandidate)
    #num_units=len(nodes)
    total_pop=sum(x[3] for x in nodes)
    #print total_pop,nodes[:10]
	#sum(nodes[x][3] for x in range(num_units))
    node_groups=[-1 for x in range(num_units)]
    if location_problem==0:
        fixed_cost_obj=0
        max_num_facility=num_districts
    if fixed_cost_obj==0:
        facilityCost=[0 for x in range(num_districts)]
    if location_problem==1 and max_num_facility<1:
        #max_num_facility=num_districts
        adaptive_number_of_facilities=1
    if location_problem==2:
        if all_units_as_candadate_locations==1:
            facilityCandidate=[x for x in range(num_districts)]
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
        if all_units_as_candadate_locations==0:
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
    if location_problem==3: #pmp
        #num_districts=num_units
        #facilityCandidate=[x for x in range(num_districts)]
        facilityCost=[0.0 for x in range(num_districts)]
        facilityCapacity=[total_pop for x in range(num_districts)]
    centersID=facilityCandidate[:]
    num_facilityCandidate=len(facilityCandidate)
    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    total_supply=sum(facilityCapacity)
    #arcpy_print("total demand: "+str(total_pop))
    #arcpy_print("total supply: "+str(total_supply))
    #arcpy_print("avg. distance to nearest facility: "+str(avg_dis_min))

    objective_overload=MAXNUMBER
    obj_balance=MAXNUMBER
    objective=MAXNUMBER
    biobjective=MAXNUMBER
    all_solutions=[]

    #best solution in each start
    best_solution =[] # node_groups[:]
    best_centersID=[]
    best_biobjective=MAXNUMBER
    best_objective=MAXNUMBER
    best_objective_overload = MAXNUMBER

    #global best solution 
    best_solution_global=[]
    best_centersID_global=[]
    best_biobjective_global = MAXNUMBER
    best_objective_global = MAXNUMBER
    best_overload_global = MAXNUMBER
    avg_dis_min =sum(nodedik[x][0] for x in range(num_units))/total_pop
    find_NearFacilityList(num_districts)
    if linear_relaxation==1:
        lplocs,sol=location_model_linear_relexation(max_num_facility,0,heuristic_time_limit,0.0001)	
        potential_facilities=[x for x in range(num_districts) if lplocs[x]>0.0001]
        print "Potential facilities by Linear Relax",potential_facilities    
    max_exclusion_list=[0.0 for x in range(num_districts)]
    
def find_NearFacilityList(nnn):
    global NearFacilityList
    if len(NearFacilityList)>0: return
    NearFacilityList=[]
    n=nnn#num_districts
    if n>num_districts: n=num_districts
    dis=0.0
    print "find_NearFacilityList()",
    for i in range(num_units):
        if i%100==0: print ".",
        fdlist=[ [x,nodedij[i][x]] for x in range(num_districts)]
        fdlist.sort(key=lambda x:x[1])
        flist=[x[0] for x in fdlist[0:n]]
        NearFacilityList.append(flist[:])
            
def cflp_mst():
    vars=[]
    if len(node_groups)<num_units: return vars
    for i in range(num_units):
        k=node_groups[i]
        if k<0: continue
        v='x_'+str(i)+ '_'+ str(k)
        vars.append([v,1])
    for i in range(num_districts):
        if centersID[i]<0:continue
        v='y_'+str(i)
        vars.append([v,1])
    return vars

def get_cand_locations(dlist,ulist):
    clist=[]
    if len(max_exclusion_list)<=0: 
        clist=[x for x in range(num_districts) if centersID[x]<0 and nearCustomer[x] in ulist]
        random.shuffle(clist)
    else:
        clist=[[x,max_exclusion_list[x]] for x in range(num_districts) if centersID[x]<0 and nearCustomer[x] in ulist and x in potential_facilities]
        clist.sort(key=lambda x:x[1])
        clist=[x[0] for x in clist]
    nc=len(clist)
    nf=len(dlist)
    nu=len(ulist)
    t=time.time()
    mnsize=min(3000,num_units*num_districts/10)
    ns=min(nc,mnsize/nu-nf)
    #ns=nc
    if ns<nf: ns=min(nc,nf)
    if ns==nc: return [clist]
    #cloclist=clist[:ns/2]
    #while len(cloclist)<ns:
    #    r=random.random()
    #    idx=int(r**1.5*0.999*(nc-ns/2))
    #    k=clist[ns/2+idx]
    #    if k not in cloclist: cloclist.append(k)
    #cloclist.append(clist[:ns])
    cloclist=[]
    clist1=clist[:ns/2]
    clist2=clist[ns/2:]
    for i in range(5):
        if initial_solution_method!=9:
            random.shuffle(clist)
            cloclist.append(clist[:ns])
        else:
            random.shuffle(clist2)
            cloclist.append(clist1+clist2[:ns/2])
    #print nf,nu,mnsize,nc+nf,ns+nf,
    return cloclist

def loc_sub_mst(dlist,ulist):
    vars=[]
    for i in ulist:
        k=node_groups[i]
        if k<0: continue
        v='x_'+str(i)+ '_'+ str(k)
        vars.append([v,1])
    for i in dlist:
        if centersID[i]<0:continue
        v='y_'+str(i)
        vars.append([v,1])
    return vars
#def pmp_mip_model(max_radius,cover_pecentage,time_limit,mipgap): 
def pmp_mip_model(num_f,radius,cover_percent,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pmp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    ulist=range(num_units)
    xvariables={}
    xcost={}
    for i in ulist:
        for k in centers:
            xvariables["x_" +str(i)+"_"+str(k)]=pulp.LpVariable("x_" +str(i)+"_"+str(k), 0, 1, pulp.LpBinary)
            xcost["x_" +str(i)+"_"+str(k)]=nodedik[i][k]*envy_objective_weight_of_travelcost
            d=nodedij[i][k]
            if envy_service_objective==1 and d>envy_service_distance:
                xcost["x_" +str(i)+ "_"+ str(k)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            if envy_service_objective==4 and d>envy_service_distance:
                xcost["x_" +str(i)+ "_"+ str(k)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
    costs={}
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        #costs["y_" +str(i)]=100000000#facilityCost[i]
    fvariables={}

    #zvariables={}    
    #for i in ulist:
    #    zvariables[i]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpBinary)

    obj=""
    for x in xvariables:
        obj +=xcost[x]* xvariables[x]
    prob += obj

    for i in ulist:
        s=""
        for k in centers:
            s+=xvariables["x_" +str(i)+"_"+str(k)]
        prob +=s  == 1
    s=""
    for k in centers:
        s+=yvariables["y_"+str(k)]
    prob +=s  == num_f

    for k in centers:
        s=""
        for i in ulist:
            s+=xvariables["x_" +str(i)+"_"+str(k)]
        s-=num_units*yvariables["y_" + str(k)]
        prob +=s <= 0

    for k in centers:
        for i in ulist:
            s=xvariables["x_" +str(i)+"_"+str(k)] -yvariables["y_" + str(k)]
            prob +=s <= 0

    prob.writeLP("_pmp.lp")
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        ##noprint "model unsolved..."
        return []
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
            if items[0]=='x':
                k=int(items[2])
                node_groups[i]=k
    '''for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k] >=0:
                node_groups[i]=k
                break'''
    #update_district_info()
    update_envy_district_info()
    return 1

def cpmp_mip_model(num_f,radius,cover_percent,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pmp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    ulist=range(num_units)
    xvariables={}
    xcost={}
    for i in ulist:
        for k in centers:
            xvariables["x_" +str(i)+"_"+str(k)]=pulp.LpVariable("x_" +str(i)+"_"+str(k), 0, 1, pulp.LpBinary)
            xcost["x_" +str(i)+"_"+str(k)]=nodedik[i][k]*envy_objective_weight_of_travelcost
            d=nodedij[i][k]
            if envy_service_objective==1 and d>envy_service_distance:
                xcost["x_" +str(i)+ "_"+ str(k)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            if envy_service_objective==4 and d>envy_service_distance:
                xcost["x_" +str(i)+ "_"+ str(k)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
    fvariables={}
    if spatial_contiguity==1:
        for i in range(num_units):
            for j in node_neighbors[i]:
                for k in range(num_districts):
                    fvariables["f_" +str(i)+ "_"+str(j)+ "_"+ str(k)]=pulp.LpVariable("f_" +str(i)+ "_"+ str(j)+ "_"+ str(k), 0, None, pulp.LpInteger)
    #zvariables={}    
    #for i in ulist:
    #    zvariables[i]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpBinary)

    obj=""
    for x in xvariables:
        obj +=xcost[x]* xvariables[x]
    prob += obj

    for i in ulist:
        s=""
        for k in centers:
            s+=xvariables["x_" +str(i)+"_"+str(k)]
        prob +=s  == 1
    s=""
    for k in centers:
        s+=yvariables["y_"+str(k)]
    prob +=s  == num_f

    for k in centers:
        s=""
        for i in ulist:
            s+=nodes[i][3]*xvariables["x_" +str(i)+"_"+str(k)]
        s-=facilityCapacity[k]*yvariables["y_" + str(k)]
        prob +=s <= 0

    prob.writeLP("_cpmpr.lp")
    #initvalues=pmp_mst(dlist,ulist)
    #for x,v in initvalues:
    #    if x.find("x")==0:  xvariables[x].setInitialValue(v)
    #    if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    solver_message=1
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        ##noprint "model unsolved..."
        return []
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
            if items[0]=='x':
                node_groups[i]=int(items[2])
    update_envy_district_info()
    return 1

def pmp_mst():
    varibles=[]
    if len(node_groups)<num_units: return varibles
    if node_groups[0]==-1:return varibles
    for i in range(num_units):
        k=node_groups[i]
        if k<0: continue
        v='x_'+str(i)+ '_'+ str(k)
        varibles.append([v,1])
    for i in range(num_districts):
        if centersID[i]<0:continue
        v='y_'+str(i)
        varibles.append([v,1])
    return varibles

def lscp_mip_model_init(max_radius,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    random.shuffle(centers)
    #centers=centers[:num_districts*9/10]
    ulist=range(num_units)
    xvariables={}
    costs={}
    yvariables={}
    
    for i in centers:
        dis=0.0
        for j in nearCustomers[i]:
            if nodedij[j][i]>max_radius:
                break
            dis+=nodedik[j][i]
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=facilityCost[i] + dis*(9.5+random.random())/10
    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpInteger)
    
    obj=""
    for x in yvariables:
        obj +=costs[x]*yvariables[x]
    prob += obj

    for k in centers:
        if k in facility_inclusion_list:
            prob += yvariables["y_" +str(k)] == 1

    #cons 2
    for i in ulist:
        s=""
        for j in centers:
            if nodedij[i][j]>max_radius: continue
            s+=yvariables["y_" + str(j)]
        prob +=s - zvariables["z_" +str(i)] >= 0

#    for i in ulist:
#        for j in centers:
#            if nodedij[i][j]>max_radius: continue
#            s=yvariables["y_" + str(j)]+ zvariables["z_" +str(i)]
#            prob +=s <= 1

    s=""
    for i in ulist:
        s+=nodes[i][3]*zvariables["z_" +str(i)]
    prob += s>=total_pop*cover_pecentage/100
    #maxuc=total_pop-total_pop*cover_pecentage/100
    #prob += s<=maxuc
    prob.writeLP("_lscp.lp")
    #initvalues=pmp_mst(dlist,ulist)
    #for x,v in initvalues:
    #    if x.find("x")==0:  xvariables[x].setInitialValue(v)
    #    if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    #solver_message=1
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        print "model unsolved..."
        return 0
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k] >=0:
                node_groups[i]=k
                break
    update_district_info()
    return 1

#TESTING
def clscp_mip_model_init(max_radius,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    global district_info
    prob = pulp.LpProblem("clscp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    sampling=0
    if total_supply>10*total_pop:
        sampling=1
        min_nf=int(total_pop*num_districts/total_supply)
        clist=range(num_districts)
        random.shuffle(clist)
        clist=clist[:min_nf*10]
        centers=list(set(clist+facility_inclusion_list))

    covered_units=[[] for x in range(num_districts)]
    covered_cost=[facilityCost[x] for x in range(num_districts)]
    for k in range(num_districts):
        covered=[]
        if k not in centers: continue
        supply=facilityCapacity[k] #*1.15 
        if sampling==0:
            supply=facilityCapacity[k]*(1+random.random()/20) #(1.15+random.random()/20)  #pie*r*r / 2.6*r*r = pie/2.6=
        cost=0.0
        for i in nearCustomers[k]:
            if nodedij[i][k]>max_radius: break
            if supply<nodes[i][3]: break
            supply-=nodes[i][3]
            covered.append(i)
            cost+=nodedik[i][k]
        covered_units[k]=covered[:]
        covered_cost[k]+=cost
    avg_covered_cost=sum(covered_cost)/num_districts
    ulist=range(num_units)
    yvariables={}
    cost={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        cost["y_" +str(i)]=facilityCost[i]
    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpInteger)
    obj=""
    #for x in yvariables:
    #    obj += yvariables[x]
    for k in centers:
        obj += covered_cost[k]*yvariables["y_" +str(k)]
        #obj += (100000000 + abs(covered_cost[k]-avg_covered_cost))*yvariables["y_" +str(k)]
    prob += obj
    s=""
    for k in centers:
        s+=facilityCapacity[k]*yvariables["y_" +str(k)]
    prob +=s >= total_pop

    for k in centers:
        if k in facility_inclusion_list:
            prob += yvariables["y_" +str(k)] == 1

    if adaptive_number_of_facilities==0:
        s=""
        for k in centers:
            s+=yvariables["y_" +str(k)]
        prob +=s == max_num_facility
    for i in ulist:
        s=""
        for j in centers:
            if i in covered_units[j]: 
                s+=yvariables["y_" + str(j)]
        s -= zvariables["z_" +str(i)]
        prob += s>= 0

    if cover_pecentage >=99: 
        for i in ulist:
            prob += zvariables["z_" +str(i)] == 1
    else:
        s=""
        for i in ulist:
            s+=nodes[i][3]*zvariables["z_" +str(i)]
        prob += s>=total_pop*cover_pecentage/100.0
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved..."
        return 0
    centersID=[-1 for x in range(num_districts)]
    node_groups=[-1 for x in range(num_units)]
    fselected=[]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
                fselected.append(i)
    for k in fselected:
        for i in covered_units[k]:
            if node_groups[i]==-1: node_groups[i]=k
            if node_groups[i]>=0: 
                k2=node_groups[i]
                if nodedij[i][k] <nodedij[i][k2]: node_groups[i]=k
    update_district_info()
    ulist=[x for x in range(num_units) if node_groups[x]==-1]
    random.shuffle(ulist)
    for i in ulist:
        for k in NearFacilityList[i]:
            if centersID[k]<0: continue
            if facilityCapacity[k]-district_info[k][1]>=nodes[i][3]:
                node_groups[i]=k
                district_info[k][1] +=nodes[i][3]
                #print ".",
                break
    ulist=[x for x in range(num_units) if node_groups[x]==-1]
    random.shuffle(ulist)
    for i in ulist:
        for k in NearFacilityList[i]:
            if centersID[k]<0: continue
            node_groups[i]=k
            #print "&",
            break
    update_district_info()
    #for x in district_info: 
    #    if x[0]>0: print x
    return 1
    
def lscp_mip_model_init2(maxr2,maxr,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    global district_info
    prob = pulp.LpProblem("_clscp",pulp.LpMinimize)
    centers=range(num_districts)
    random.shuffle(centers)
    #centers=centers[:num_districts*9/10]
    centers=list(set(centers+facility_inclusion_list))
    ulist=range(num_units)
    yvariables={}
   
    costs={}
    for i in centers:
        dis=0.0
        for j in nearCustomers[i]:
            if nodedij[j][i]>maxr2:
                break
            dis+=nodedik[j][i]
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=facilityCost[i] + dis*(19.5+random.random())/20

    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpInteger)
    obj=""
    for x in yvariables:
        obj += costs[x]*yvariables[x]
    prob += obj

    for k in centers:
        if k in facility_inclusion_list:
            prob += yvariables["y_" +str(k)] == 1

    for i in ulist:
        s=""
        for j in centers:
            if nodedij[i][j]>maxr: continue
            s+=yvariables["y_" + str(j)]
        s -= zvariables["z_" +str(i)]
        prob += s>= 0

    for i in ulist:
        s=""
        for j in centers:
            if nodedij[i][j]>maxr2: continue
            s+=yvariables["y_" + str(j)]
        prob += s >= 1

    s=""
    for i in ulist:
        s+=nodes[i][3]*zvariables["z_" +str(i)]
    prob += s>=total_pop*cover_pecentage/100.0
    
    prob.writeLP("_sclp_init3.lp")
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved..."
        return 0
    centersID=[-1 for x in range(num_districts)]
    node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
    for i in ulist:
        for k in NearFacilityList[i]:
            if centersID[k]<0: continue
            node_groups[i]=k
            break
    return 1   
 
def clscp_mip_model_init2(maxr2,maxr,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    global district_info
    prob = pulp.LpProblem("_clscp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    sampling=0
    if total_supply>10*total_pop:
        sampling=1
        min_nf=int(total_pop*num_districts/total_supply)
        clist=range(num_districts)
        random.shuffle(clist)
        clist=clist[:min_nf*10]
        centers=list(set(clist+facility_inclusion_list))
    covered_units=[[] for x in range(num_districts)]
    covered_cost=[facilityCost[x] for x in range(num_districts)]
    for k in range(num_districts):
        covered=[]
        if k not in centers: continue
        supply=facilityCapacity[k]*1.15  #pie*r*r / 2.6*r*r = pie/2.6=
        #if sampling==0:
        #    supply=facilityCapacity[k] *(1+random.random()/20)  #pie*r*r / 2.6*r*r = pie/2.6=
        cost=0.0
        for i in nearCustomers[k]:
            if supply<nodes[i][3]: continue
            supply-=nodes[i][3]
            covered.append(i)
            cost+=nodedik[i][k]
            if nodedij[i][k]>maxr2:
                cost+=nodes[i][3]*penalty_on_demand
            if envy_service_objective==1 and nodedij[i][k]>envy_service_distance:
                cost+=nodes[i][3]*envy_objective_weight*(nodedij[i][k]-envy_service_distance)*(nodedij[i][k]-envy_service_distance)
        covered_units[k]=covered[:]
        covered_cost[k]+=cost
    avg_covered_cost=sum(covered_cost)/num_districts
    ulist=range(num_units)
    yvariables={}
    #for k in range(num_districts):
    #    print k, facilityCapacity[k],len(covered_units[k])
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpBinary) #maxr2 covered
    for i in ulist:
        zvariables["z2_" +str(i)]=pulp.LpVariable("z2_" +str(i), 0, 1, pulp.LpBinary) ##maxr covered

    zvariable=pulp.LpVariable("z", 0, None, pulp.LpInteger)
    obj=""
    for k in centers:
        obj += covered_cost[k]*yvariables["y_" +str(k)] 
    for i in ulist:
        obj += nodes[i][3]*penalty_on_demand*zvariables["z2_" +str(i)] 
    obj += penalty_on_demand * zvariable
    prob += obj
    '''
    s=""
    for k in centers:
        s+=facilityCapacity[k]*yvariables["y_" +str(k)]
    prob +=s >= total_pop'''

    #for k in centers:
    #    if k in facility_inclusion_list:
    #        prob += yvariables["y_" +str(k)] == 1

    if adaptive_number_of_facilities==0:
        s=""
        for k in centers:
            s+=yvariables["y_" +str(k)]
        prob +=s == max_num_facility

    for i in ulist:
        s=""
        for j in centers:
            if i in covered_units[j]: 
                s+=yvariables["y_" + str(j)]
        s+=zvariables["z2_" +str(i)]
        prob += s>= 1

    for i in ulist:
        s=""
        for j in centers:
            if i in covered_units[j]: 
                if nodedij[i][j]<=maxr: 
                    s+=yvariables["y_" + str(j)]
        s -= zvariables["z_" +str(i)]
        prob += s>= 0

    s=""
    for i in ulist:
        s+=nodes[i][3]*zvariables["z_" +str(i)]
    s+=zvariable
    prob += s>=total_pop*cover_pecentage/100.0

    prob.writeLP("_sclp_init3.lp")
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved..."
        return 0
    centersID=[-1 for x in range(num_districts)]
    node_groups=[-1 for x in range(num_units)]
    fselected=[]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            if items[0]=='z': continue
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
                fselected.append(i)
    for k in fselected:
        for i in covered_units[k]:
            if node_groups[i]==-1: node_groups[i]=k
            if node_groups[i]>=0: 
                k2=node_groups[i]
                if nodedij[i][k] <nodedij[i][k2]: node_groups[i]=k
    update_district_info()
    ulist=[x for x in range(num_units) if node_groups[x]==-1]
    random.shuffle(ulist)
    for i in ulist:
        for k in NearFacilityList[i]:
            if centersID[k]<0: continue
            if facilityCapacity[k]-district_info[k][1]>=nodes[i][3]:
                node_groups[i]=k
                district_info[k][1] +=nodes[i][3]
                break
    ulist=[x for x in range(num_units) if node_groups[x]==-1]
    random.shuffle(ulist)
    for i in ulist:
        for k in NearFacilityList[i]:
            if centersID[k]<0: continue
            node_groups[i]=k
            #print "&",
            break
    update_district_info()
    return 1

def update_cflpr_district_info(max_radius2,max_radius,cover_pecentage):
    global biobjective
    global objective_pentalty_on_rmax
    global objective_pentalty_on_covering  
    global objective_pentalty_on_overload
    global objective_envy
    global objective_covered
    global objective_rmax_not_covered
    update_district_info()
    objective_pentalty_on_overload=objective_overload*penalty_on_demand
    covered=0
    covered2=0
    for i in range(num_units):
        k=node_groups[i]
        if nodedij[i][k]<=max_radius: covered+=nodes[i][3]
        if nodedij[i][k]>max_radius2: covered2+=nodes[i][3]
    objective_covered=covered
    objective_rmax_not_covered=covered2
    objective_pentalty_on_covering= max(0, total_pop*cover_pecentage/100-covered)*penalty_on_demand
    objective_pentalty_on_rmax=covered2*penalty_on_demand
    
    #penalty_demands =objective_overload
    #penalty_demands+= max(0, total_pop*cover_pecentage/100-covered)
    #penalty_demands+=covered2
    #penalty=penalty_on_demand
    #biobjective=objective+penalty*penalty_demands
    biobjective=objective*envy_objective_weight_of_travelcost+objective_pentalty_on_overload+objective_pentalty_on_rmax+objective_pentalty_on_covering
    if fixed_cost_obj==1: biobjective += objective_fcost
    if envy_service_objective==1:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
        objective_envy=envyobj
        biobjective += envyobj
    if envy_service_objective==2:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*(d-envy_service_distance)*(d-envy_service_distance)
        var=envy_coefficient_variation*envy_service_distance
        maxobj= 0.5* total_pop * var*var
        objective_envy=envyobj
        if envyobj>maxobj: 
            biobjective +=penalty_on_demand*100*max(0,envyobj-maxobj)
    if envy_service_objective==4:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
        objective_envy=envyobj
        biobjective += envyobj
        
def update_envy_district_info():
    global biobjective
    global objective_envy
    update_district_info()
    penalty=penalty_on_demand
    if envy_service_objective==1:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
        objective_envy=envyobj
        biobjective += envyobj
    if envy_service_objective==4:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
        objective_envy=envyobj
        biobjective += envyobj
    if envy_service_objective==2:
        envyobj=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envyobj+=nodes[i][3]*(d-envy_service_distance)*(d-envy_service_distance)
        var=envy_coefficient_variation*envy_service_distance
        maxobj= 0.5* total_pop * var*var
        objective_envy=0.0
        if envyobj>maxobj: 
            objective_envy=penalty*100*(envyobj-maxobj)
        biobjective += objective_envy

def update_cflpr_centers(r,p):
    global node_groups
    global centersID
    global time_update_centers
    t=time.time()
    centers=[x for x in range(num_districts) if centersID[x]>=0 and x not in facility_inclusion_list]
    random.shuffle(centers)
    for k in centers:
        newk=update_cflpr_center(k,r,p)
        if newk==k: continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        for x in ulist: node_groups[x]=newk
        centersID[k]=-1
        centersID[newk]=facilityCandidate[newk]
        #print "a center updated!"
    time_update_centers+=time.time()-t


def update_cflpr_center(k,r,p):
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    covered=0
    for i in ulist:
        j=node_groups[i]
        if nodedij[i][k]<r: 
            covered+=nodes[i][3]
    demand=sum(nodes[x][3] for x in ulist)
    best_cost=sum(nodedik[x][k] for x in ulist)
    best_center=k
    uid=nearCustomer[k]
    centers=[]
    for x in NearFacilityList[uid]:
        if centersID[x]<0: centers.append(x)
        if len(centers)>=min(20,len(ulist)): break
    for i in centers:
        if i==k: continue
        if facilityCapacity[i]<demand and facilityCapacity[k]<facilityCapacity[i]: continue
        cost=sum(nodedik[x][i] for x in ulist)
        if cost>=best_cost: continue
        covered3=0
        for j in ulist:
            if nodedij[j][i]<r:covered3+=nodes[j][3]
        if covered3<covered: continue 
        best_cost=cost
        best_center=i
    return best_center


def update_zero_demand():
    global node_groups
    for i in range(num_units):
        if nodes[i][3]>0: continue
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break

##LBc << LBr, constrained by radius, using SCLP to generate solutions.
##LBc >> LBr, constrained by capacity, using cSCLP to generate solutions.
##LBc == LBr, constrained by both radius and capacity? using cSCLP to generate solutions? drop???
def cflpr_matheuristic(max_radius,radius_preference,cover_pecentage,multi_start,maxloops):
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global best_centersID_global
    global best_solution_global
    global objective
    global biobjective
    global objective_overload
    global node_groups
    global centersID
    global district_info
    global facilityCost
    #global node_neighbors
    global region_pool
    global pool_index
    global is_spp_modeling
    global all_solutions
    global max_loops_solution_not_improved
    global multi_start_count
    global solver_message
    global pop_dis_coeff
    global cflp_max_service_radius
    global cflp_service_radius_preference
    global cflp_min_service_coverage_percentage
    global penalty_on_demand

    initialize_instance()
    max_loops_solution_not_improved=maxloops
    multi_start_count=multi_start
    print "num_facility:",max_num_facility
    pop_dis_coeff=100000000*1000.0/total_pop
    covered_demand_nodes=[]
    opt_radius=[0.0 for x in range(num_districts)]
    cflp_max_service_radius=max_radius
    cflp_service_radius_preference=radius_preference
    cflp_min_service_coverage_percentage=cover_pecentage
    #penalty_on_demand=10000*sum(facilityCost)/sum(facilityCapacity)
    penalty_on_demand=current_penalty()
    if fixed_cost_obj==0.0: 
        penalty_on_demand=10000.0
    if envy_service_objective==1: penalty_on_demand*=envy_objective_weight
    if envy_service_objective==4: penalty_on_demand*=envy_objective_weight
    print "penalty_on_demand:",penalty_on_demand
    all_solutions=[]
    region_pool=[]
    t=time.time()
    best_biobjective_global = MAXNUMBER
    best_biobjective = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[[] for x in range(num_districts)]
    if sum(facilityCost)<1:
        for x in range(num_districts):
            facilityCost[x]=100000000
    if fixed_cost_obj==0:
        for x in range(num_districts):
            facilityCost[x]=0
    penalty_on_demand=10000.0

    not_improved_global=0
    init_methods=[0,1,2,3,4]
    #if adaptive_number_of_facilities==0: init_methods=[1,2]
    loops=multi_start_count
    method=initial_solution_method
    if method<0: method=0
    for idx in range(loops*2):
        init_t=time.time()
        sta=0
        if method==0: #cflp
            drop_method_cflpr(max_radius,radius_preference,cover_pecentage)
        if method==1: #DC-cflp
            sta=clscp_mip_model_init2(max_radius,radius_preference,cover_pecentage,100,0.01)
        if method==2: #c-cflp
            sta=clscp_mip_model_init(radius_preference,cover_pecentage,100,0.01)
        if method==3: #LSCP with Rmax
            sta=1
            clist=range(num_districts)
            for x in clist: centersID[x]=-1
            random.shuffle(clist)
            n=int(total_pop*num_districts/sum(facilityCapacity))+1
            for x in clist[:n]:
                centersID[x]=facilityCandidate[x]
            for i in range(num_units):
                for k in NearFacilityList[i]:
                    if centersID[k]>=0:
                        node_groups[i]=k
                        break
        if sta<1: 
            method=0
            sta=drop_method_cflpr(max_radius,radius_preference,cover_pecentage)
        if penalty_on_demand==0.0:
            update_district_info()
            penalty_on_demand=10000*objective/total_pop
            print "new penalty_on_demand:",penalty_on_demand
        update_zero_demand()
        #check_required_facility()
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
        '''assign_CFLPr1(max_radius,radius_preference,cover_pecentage)
        assign_CFLPr2(max_radius,radius_preference,cover_pecentage)
        assign_CFLPr3(max_radius,radius_preference,cover_pecentage)'''
        VND_local_search()
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)

        coverdemand=0
        coverdemand2=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=radius_preference: coverdemand+=nodes[i][3]
            if nodedij[i][k]>max_radius: coverdemand2+=nodes[i][3]
        cp=coverdemand*1.0/total_pop
        penalty_demand=objective_overload+coverdemand2+max(total_pop*cover_pecentage/100-coverdemand,0)
        #penalty=100*sum(facilityCost)/sum(facilityCapacity)
        penalty=penalty_on_demand
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
        nf=sum(1 for x in centersID if x>=0)
        update_region_pool_all()
        update_best_solution()
        print "init. sol.",idx, "mthd",method,nf,biobjective,objective,objective_overload,coverdemand2, int(cp*100)/100.0,time.time()-t
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        sol=print_equality_measures()
        s="stat: "+str(sol["Mean"]) +" "+str(sol["MeanDev"]) +" "+str(sol["StdDev"]) +" "+str(sol["Gini"])
        print s
        if len(population)>=multi_start_count: break
    population.sort(key=lambda x:x[0])
    all_solutions=population
    print "pop objs: "+str(population[0][0])+"-"+str(population[-1][0])
    loop=-1
    while 1:
        r=random.random()
        sidx = int(min(multi_start_count-1,len(population))* pow(r,1.5)*0.999)
        #if not_improved_global>= max_loops_solution_not_improved/5:
        #    sidx=0
        node_groups=population[sidx][2][:]
        centersID=population[sidx][1][:]
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
        bnf=sum(1 for x in centersID if x>=0)
        orgs=str(bnf) +" "+ str(int(biobjective))+" " +str(int(objective))+" "+str(objective_overload)+" -> "
        obj_old=biobjective
        loop+=1
        not_improved_global+=1
        obj=best_biobjective_global
        sscflpr_sub_model(max_radius,radius_preference,cover_pecentage,20,0.000000001)
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
        update_best_solution()
        update_region_pool_all()

        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        population.sort(key=lambda x:x[0])
        population=pop_selection(population)
        all_solutions=population

        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=radius_preference: coverdemand+=nodes[i][3]
        s="" 
        if best_biobjective_global<obj-0.000001:  #0.1%
            s="*"
            impp=(best_biobjective_global-obj)/(best_biobjective_global-best_objective_fcost_global)*1000 #-0.1%
            not_improved_global+=int(max_loops_solution_not_improved*impp)
            if not_improved_global<0: not_improved_global=0
        else: s="-"
        bnf=sum(1 for x in best_centersID_global if x>=0)
        s+="Loop "+str(loop) + " Best: " +str(bnf)+" "+str(int(best_biobjective_global))+" "+str(int(best_objective_global))+" "+str(int(best_overload_global))
        bnf=sum(1 for x in centersID if x>=0)
        s+=" Current: "+ orgs + str(bnf)+" "+str(int(biobjective))+" "+str(int(objective))+" "+str(int(objective_fcost))
        s+=" Info: " + str(not_improved_global)+ " "+str(int(time.time()-t))
        coverdemand=0
        coverdemand2=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=radius_preference: coverdemand+=nodes[i][3]
            if nodedij[i][k]>max_radius: coverdemand2+=nodes[i][3]
        cp=coverdemand*1.0/total_pop
        penalty_demand=objective_overload+coverdemand2+max(total_pop*cover_pecentage/100-coverdemand,0)
        s+=" "+str(int(objective_overload)) + " "+  str(coverdemand2) +" "+ str(max(total_pop*cover_pecentage/100-coverdemand,0)) +" "+ str(int(cp*10000)/100.0)
        s+=" "+str(int(biobjective-obj_old)) #+" "+str(population[0][0])+"-"+str(population[-1][0])
        arcpy_print(s)
     
        if not_improved_global >= max_loops_solution_not_improved: break
    #post procedure
    ga_time=time.time()-t
    node_groups=best_solution_global[:]
    centersID=best_centersID_global[:]
    update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
    nf=sum(1 for x in centersID if x>=0)
    print "Heuristic solution:",nf,biobjective,objective_fcost,objective,objective_overload,ga_time
    t_spp=time.time()
    if is_spp_modeling>=1:
        arcpy_print("SPP modelling..."+str(len(region_pool)) )
        sta=cflpr_scp_model(max_radius,radius_preference,cover_pecentage,ga_time/10,0.00000001)
        if sta>0:
            update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
            update_best_solution()
            nf=sum(1 for x in centersID if x>=0)
            print "spp solution:",nf,biobjective,objective_fcost,objective,objective_overload,objective/total_pop,time.time()-t_spp
            population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,max(0,total_pop-objective_supply)])
    penalty_on_demand=max(max(facilityCost),10000)
    for idx in range(len(all_solutions)):
        node_groups=all_solutions[idx][2][:]
        centersID=all_solutions[idx][1][:]
        update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
        all_solutions[idx][0]=biobjective
    population.sort(key=lambda x:x[0])

    node_groups=population[0][2][:]
    centersID=population[0][1][:]
    update_cflpr_district_info(max_radius,radius_preference,cover_pecentage)
    coverdemand=0
    for i in range(num_units):
        k=node_groups[i]
        if nodedij[i][k]<=radius_preference: coverdemand+=nodes[i][3]
    nf=sum(1 for x in centersID if x>=0)
    print "final solution:",nf,biobjective,objective_fcost,objective,objective_overload,objective/total_pop,int(coverdemand*10000/total_pop)/100.0
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

def current_penalty():
    #return 10*sum(facilityCost)/sum(facilityCapacity)
    return 100*max(facilityCost)*num_districts/sum(facilityCapacity)

def sscflpr_sub_model(max_radius2,max_radius,cover_pecentage,time_limit,mipgap): 
    #return sscflpr_sub_model2(max_radius2,max_radius,cover_pecentage,time_limit,mipgap)
    #may be more facility needed, due to the covering constraint 
    global centersID
    global node_groups
    rand=0
    if random.random()>10.5: rand=1
    dlist=[]
    ulist=[]


    #dlist=select_region(u)
    #dlist=select_region_adaptive(u)
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    r=random.random()
    neig=1
    if r<0.25:
        #dlist=select_region_verylarge(u)
        dlist=select_region(-1)
        ulist=[x for x in range(num_units) if nodes[x][3]>0 and node_groups[x] in dlist]
        clist=get_cand_cflpr_locations(dlist,ulist,2)
        neig=1
    elif r<0.5:
        dlist=select_region_adaptive(-1)
        ulist=[x for x in range(num_units) if nodes[x][3]>0 and node_groups[x] in dlist]
        clist=get_cand_cflpr_locations(dlist,ulist,2)
        neig=2
    elif r<0.75:
        if nf<=30:
            dlist=[x for x in range(num_districts) if centersID[x]>=0]
            ulist=range(num_units)
        else:
            u=random.randint(0,num_units-1)
            #dlist=[x for x in NearFacilityList[u] if centersID[x]>=0]
            #dlist=dlist[:30]
            dlist=[]
            u=random.randint(0,num_units-1)
            for k in NearFacilityList[u]:
                if centersID[k]>=0:
                    dlist.append(k)
                if len(dlist)>=30:
                    break
            ulist=[x for x in range(num_units) if nodes[x][3]>0 and node_groups[x] in dlist]
        #clist=[x for x in range(num_districts) if centersID[x]<0]
        clist=[NearFacilityList[x][0] for x in ulist]
        clist=[x for x in clist if centersID[x]<0]
        clist=list(set(clist))
        random.shuffle(clist)
        clist=clist[:min(nf,30)/3] #[:10+len(dlist)/10]
        neig=3
    else:
        dlist=select_corridor_region()
        ulist=[x for x in range(num_units) if nodes[x][3]>0 and node_groups[x] in dlist]
        clist=get_cand_cflpr_locations(dlist,ulist,2)
        neig=4
    #print neig,
  
    centers=list(set(clist+dlist))
    #print [len(dlist),len(centers)],
    exist_facility=[x for x in range(num_districts) if centersID[x]>=0 and x not in centers and facilityCapacity[x]-district_info[x][1]>0]
    #exist_facility=[]
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    xvariables={}
    costs={}

    #penalty=current_penalty()
    penalty=penalty_on_demand
    #print int(penalty),
    
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j] *envy_objective_weight_of_travelcost/50
            d=nodedij[i][j]
            if d>max_radius2:
                costs["x_" +str(i)+ "_"+ str(j)]+= nodes[i][3]*penalty
            if envy_service_objective==1 and d>envy_service_distance:
                costs["x_" +str(i)+ "_"+ str(j)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            if envy_service_objective==4 and d>envy_service_distance:
                costs["x_" +str(i)+ "_"+ str(j)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
            if rand==1: costs["x_" +str(i)+ "_"+ str(j)]*=(99.5+random.random())/100
        for j in exist_facility:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]*envy_objective_weight_of_travelcost/50
            d=nodedij[i][j]
            if d>max_radius2:
                costs["x_" +str(i)+ "_"+ str(j)] += nodes[i][3]*penalty
            if envy_service_objective==1 and d>envy_service_distance:
                costs["x_" +str(i)+ "_"+ str(j)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            if envy_service_objective==4 and d>envy_service_distance:
                costs["x_" +str(i)+ "_"+ str(j)]+=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
            if rand==1: costs["x_" +str(i)+ "_"+ str(j)]*=(99.5+random.random())/100
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=facilityCost[i]
        if rand==1:
            costs["y_" +str(i)]*=(29.5+random.random())/30
    zvariable=pulp.LpVariable("z", 0, None, pulp.LpContinuous) #cover gap
    if envy_service_objective==2:
        z3variable=pulp.LpVariable("z3", 0, None, pulp.LpContinuous) #envy gap
    z2variables={}
    for i in centers:
        z2variables["z2_" +str(i)]=pulp.LpVariable("z2_" +str(i), 0, None, pulp.LpInteger)

    obj=""
    for x in yvariables:
        obj +=costs[x]* yvariables[x]
    for x in xvariables:
            obj+=costs[x]*xvariables[x]
    for x in z2variables:
        obj += penalty*z2variables[x]
    obj+=penalty*zvariable
    if envy_service_objective==2:
        obj+= penalty*10000*z3variable
    prob += obj

    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        for j in exist_facility:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        prob += s==1

    for k in centers:
        s=""
        for i in ulist:
            s+= nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(k) ]
        s-= facilityCapacity[k]*yvariables["y_" +str(k)]
        s-= z2variables["z2_" +str(k)]
        prob += s <= 0

    for k in exist_facility:
        s=""
        for i in ulist:
            s+=nodes[i][3]*xvariables["x_" +str(i)+ "_"+ str(k)]
        prob+= s <= facilityCapacity[k] - district_info[k][1]

    for k in centers:
        s= z2variables["z2_" +str(k)]
        s-= facilityCapacity[k]*yvariables["y_" +str(k)]
        prob += s <= 0

    for k in centers:
        if k in facility_inclusion_list and envy_service_objective==0:
            prob += yvariables["y_" +str(k)] == 1

    total_cover1=0
    total_cover2=0
    for i in range(num_units):
        if nodedij[i][node_groups[i]]>max_radius: continue
        if i in ulist: 
            total_cover1+=nodes[i][3]
        else: total_cover2+=nodes[i][3]
    min_covered=total_pop*cover_pecentage/100-total_cover2
    s=""
    for i in ulist:
        for j in centers:
            if nodedij[i][j]>max_radius: continue
            s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
        for j in exist_facility:
            if nodedij[i][j]>max_radius: continue
            s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
    s+=zvariable
    prob += s>=min_covered

    s=""
    for x in yvariables: s+=yvariables[x]
    numf=sum(1 for x in centersID if x>=0)
    if adaptive_number_of_facilities==0:
        if numf==max_num_facility:
            prob += s== len(dlist)
    else:
        prob += s <= len(dlist)+1
        '''
        soft_demand=objective_overload
        soft_demand+= sum([nodes[x][3] for x in range(num_units) if nodedij[x][node_groups[x]] >max_radius2])
        covered=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<max_radius:
                covered+=nodes[i][3]
        soft_demand+= max(0, total_pop*cover_pecentage/100-covered)
        if soft_demand==0: 
             prob += s >= len(dlist)-1
             prob += s <= len(dlist)+1
        elif soft_demand<=total_pop/nf/100: 
             prob += s >= len(dlist)-1
             prob += s <= len(dlist)+1
        elif soft_demand<=total_pop/nf: 
             prob += s >= len(dlist)
             prob += s <= len(dlist)+2
        else:
             prob += s <= len(dlist)+2
        '''
    if envy_service_objective==2:
        envy=0.0
        for i in range (num_units):
            if i in ulist: continue
            k=node_groups[i]
            d=nodedij[i][k]
            if d<= envy_service_distance: continue
            envy+= nodes[i][3]*(d-envy_service_distance)*(d-envy_service_distance)
        var=envy_coefficient_variation*envy_service_distance
        maxobj= 0.5* total_pop * var*var
        s=""
        for i in ulist:
            for k in centers:
                d=nodedij[i][k]
                if d<= envy_service_distance: continue
                s+= nodes[i][3] * (d - envy_service_distance) *(d - envy_service_distance) * xvariables["x_" +str(i)+ "_"+ str(k)]
            for k in exist_facility:
                d=nodedij[i][k]
                if d<= envy_service_distance: continue
                s+= nodes[i][3] * (d - envy_service_distance) *(d - envy_service_distance) * xvariables["x_" +str(i)+ "_"+ str(k)]
        if len(s)>1:
            prob +=s - z3variable <= maxobj-envy
    
    prob.writeLP("_SSCFLP.lp")
    initvalues=loc_sub_mst(dlist,ulist)
    for x,v in initvalues:
        if x.find("x")==0:  xvariables[x].setInitialValue(v)
        if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap),'set mip tolerances absmipgap 100', 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit,options=[("MIPGap",gap),("MIPGapAbs",gap*objective),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    
    if prob.status<=0:
        print "model unsolved... or infeasible, msg=", prob.status
        prob.writeLP("_sscflpr_sub_model.lp")
        return 0
    
    nlist=[]
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    old_centers = centersID[:] #[x for x in range(num_districts) if centersID[x]>=0]
    for x in dlist: centersID[x]=-1
    #for x in ulist: node_groups[x]=-1
    for v in prob.variables():
        if (v.varValue >= 0.9):
            if v.name=="z": continue
            if v.name=="z3": continue
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='z2': continue
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
                nlist.append(i)
            if items[0]=='x':
                node_groups[i]=int(items[2])
    #print [neig,len(dlist),len(centers),len(nlist),old_centers!=centersID],
    '''nf2=sum(1 for x in range(num_districts) if centersID[x]>=0)
    new_centers=[x for x in range(num_districts) if centersID[x]>=0]
    if len(dlist)>len(nlist) and nf2>nf:
        dlist.sort()
        nlist.sort()
        old_centers.sort()
        new_centers.sort()
        print dlist,nlist
        print old_centers
        print new_centers'''
    return 1

def pmp_model(p,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pmp",pulp.LpMinimize)
    centers=facilityCandidate
    xvariables={}
    costs={}
    ulist=range(num_units)
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=facilityCost[i]
    
    obj=""
    for x in xvariables:
        obj += costs[x]*xvariables[x]
    prob += obj

    #con 1
    s=""
    for k in centers:
        s+=yvariables["y_" +str(k)]
    prob +=s == p

    #cons 2
    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        prob +=s == 1

    for i in ulist:
        for k in centers:
            s= xvariables["x_" +str(i) + "_"+ str(k) ] - yvariables["y_" +str(k)]
            prob +=s <= 0

    #prob.writeLP("_sub_pmp.lp")
    #initvalues=pmp_mst(dlist,ulist)
    #for x,v in initvalues:
    #    if x.find("x")==0:  xvariables[x].setInitialValue(v)
    #    if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        ##noprint "model unsolved..."
        return []
    nlist=[]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                nlist.append(i)
    centersID=[-1 for x in range(num_districts)]
    for x in nlist: centersID[x]=x
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break
    update_district_info()
    return 1

def pmp_TB():
    global node_groups
    global centersID
    best_exch=[] #klist,save
    best_obj=biobjective
    clist=[x for x in range(num_districts) if centersID[x]>=0]
    random.shuffle(clist)
    klist=[x for x in range(num_districts) if centersID[x]<0]
    random.shuffle(klist)
    for k in clist:
        #u=nearCustomer[k]
        #klist=NearFacilityList[u][num_districts/20:num_districts/5]
        random.shuffle(klist)
        for nk in klist:
            #if nk in clist: continue
            obj=0.0
            new_centersID=centersID[:]
            new_centersID[k]=-1
            new_centersID[nk]=facilityCandidate[nk]
            for i in range(num_units):
                j=node_groups[i]
                if j!=k:
                    obj+=min(nodedik[i][j],nodedik[i][nk])
                else:
                    for j in NearFacilityList[i]:
                        if new_centersID[j]>=0:
                            obj+=nodedik[i][j]
                            break
            if envy_service_objective==1:
                for i in range(num_units):
                    for j in NearFacilityList[i]:
                        if new_centersID[j]>=0:
                            d=nodedij[i][j]
                            if d>envy_service_distance:
                                obj+= nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
                            break
            if envy_service_objective==4:
                for i in range(num_units):
                    for j in NearFacilityList[i]:
                        if new_centersID[j]>=0:
                            d=nodedij[i][j]
                            if d>envy_service_distance:
                                obj+= nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
                            break
            if obj<best_obj: 
                best_obj=obj
                best_exch.append([k,nk,obj])
                #break
            #print obj,
            if len(best_exch)>=1: break
        if len(best_exch)>=1: break
    if best_exch==[]: 
        #TB_tabu_list.append(centersID[:])
        #print "tabu len=",len(TB_tabu_list)
        return 0
    #print best_exch,best_obj,biobjective
    best_exch.sort(key=lambda x:x[2])
    k=best_exch[0][0]
    centersID[k]=-1
    nk=best_exch[0][1]
    centersID[nk]=facilityCandidate[nk]
    obj=biobjective
    for i in range(num_units):
        j=node_groups[i]
        if j!=k:
            if nodedik[i][j]>nodedik[i][nk]: node_groups[i]=nk
        else:
            for j in NearFacilityList[i]:
                if centersID[j]>=0:
                    node_groups[i]=j
                    break
    #if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    update_district_info()
    #update_centers()
    #print "tb",obj,biobjective,best_exch
    #if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    return 1

def pmp_envy_find(fi,facility2): #find fr for fi
    savings=0.0
    loss=[[x,0.0] for x in range(num_districts)]
    if fixed_cost_obj==1: 
        savings-=facilityCost[fi]
        loss=[[x,-facilityCost[x]] for x in range(num_districts)]
    w=envy_objective_weight_of_travelcost
    for i in range(num_units):
        k=node_groups[i]
        if nodedik[i][fi]<nodedik[i][k]: 
            savings+=nodedik[i][k]*w-nodedik[i][fi]*w
        else:
            loss[k][1]+=min(nodedik[i][fi],nodedik[i][facility2[i]])*w-nodedik[i][k]*w
        if envy_service_objective==1:
            envy1=0.0
            envy2=0.0
            envy3=0.0
            d=nodedij[i][k]
            if d>envy_service_distance:
                envy1=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            d=nodedij[i][fi]
            if d>envy_service_distance:
                envy2=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            d=nodedij[i][facility2[i]]
            if d>envy_service_distance:
                envy3=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            if nodedij[i][fi]<nodedij[i][k]:
                savings+=envy1-envy2
            else:
                loss[k][1]+=min(envy2,envy3)-envy1
        if envy_service_objective==4:
            envy1=0.0
            envy2=0.0
            envy3=0.0
            d=nodedij[i][k]
            if d>envy_service_distance:
                envy1=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
            d=nodedij[i][fi]
            if d>envy_service_distance:
                envy2=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
            d=nodedij[i][facility2[i]]
            if d>envy_service_distance:
                envy3=nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
            if nodedij[i][fi]<nodedij[i][k]:
                savings+=envy1-envy2
            else:
                loss[k][1]+=min(envy2,envy3)-envy1
    loss=[x for x in loss if centersID[x[0]]>=0]
    loss.sort(key=lambda x:x[1])
    fr=loss[0][0]
    profit=savings-loss[0][1]
    return fr,profit
    
def pmp_find(fi,facility2): #find fr for fi
    if envy_service_objective==1:
        return pmp_envy_find(fi,facility2)
    if envy_service_objective==4:
        return pmp_envy_find(fi,facility2)
    savings=0.0
    loss=[[x,0.0] for x in range(num_districts)]
    if fixed_cost_obj==1: 
        savings-=facilityCost[fi]
        loss=[[x,-facilityCost[x]] for x in range(num_districts)]
    for i in range(num_units):
        k=node_groups[i]
        if nodedik[i][fi]<nodedik[i][k]: 
            savings+=nodedik[i][k]-nodedik[i][fi]
        else:
            loss[k][1]+=min(nodedik[i][fi],nodedik[i][facility2[i]])-nodedik[i][k]
    loss=[x for x in loss if centersID[x[0]]>=0]
    loss.sort(key=lambda x:x[1])
    fr=loss[0][0]
    profit=savings-loss[0][1]
    return fr,profit

#On the implementation of a swap-based local search procedure for the p-median problem
def pmp_Whitaker():
    #return pmp_TB()
    global node_groups
    global centersID
    global time_Whitaker
    #global TB_tabu_list
    #if centersID in TB_tabu_list: return -1
    t=time.time()
    facility2=[0 for x in range(num_units)] #num_districts
    for i in range(num_units):
        k1=node_groups[i]
        for k2 in NearFacilityList[i]:
            if centersID[k2]>=0 and k2!=k1: 
                facility2[i]=k2
                break
    klist=[x for x in range(num_districts) if centersID[x]<0]
    random.shuffle(klist)
    improved=0
    best_swap=[-1,-1,0.0]
    swap=[]
    for nk in klist:
        k,profit=pmp_find(nk,facility2)
        if profit<=0: continue
        swap.append([k,nk,profit])
        if profit>best_swap[2]:
            best_swap[0]=k
            best_swap[1]=nk
            best_swap[2]=profit
        if len(swap)>5:
            break
    swap.sort(key=lambda x:-x[2])
    if best_swap[0]>=0:
        #k=best_swap[0]
        #nk=best_swap[1]
        idx=0
        if len(swap)>2:
            idx=random.randint(0,2)
        #idx=0
        k=swap[idx][0]
        nk=swap[idx][1]
        centersID[k]=-1
        centersID[nk]=facilityCandidate[nk]
        for i in range(num_units):
            j=node_groups[i]
            if j!=k:
                if nodedik[i][j]>nodedik[i][nk]: node_groups[i]=nk
            else:
                for j in NearFacilityList[i]:
                    if centersID[j]>=0:
                        node_groups[i]=j
                        break
        update_district_info()
        #print [x for x in centersID if x>=0]
        improved=1
    time_Whitaker+=time.time()-t
    #if improved==-1: 
    #    if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    #print len(swap),
    return improved

def spp_mst():
    vars=[]        
    num_pool=len(region_pool)
    for k in range(num_districts):
        ulist=[x for x in range(num_units) if node_groups[x]==k] 
        for i in range(num_pool):
            if ulist==region_pool[i][0] and region_pool[i][4]==k:
                vars.append([i,1])
                break
    return vars

def sppmodel(maxtime,mipgap):
    global node_groups
    global centersID
    if mip_solver not in mip_solvers:
        return 0
    if len(region_pool)<=10:
        ##noprint "no candidate district!"
        print "len(region_pool)<=10", len(region_pool)
        return 0
    alpha_coeff=avg_dis_min*pop_dis_coeff  
    prob = pulp.LpProblem("_spp",pulp.LpMinimize)
    variables=[]
    costs=[]
    ###noprint "buiding spp model..."
    #cost_stat=[]
    #cost_sum=[[0,0.0,0.0] for k in range(num_districts)]
    #for x in region_pool: print x
    for i in range(len(region_pool)):
        x=pulp.LpVariable("x_" +"{0:07d}".format(i), 0, 1,pulp.LpBinary)
        variables.append(x)
        cost=region_pool[i][1]*envy_objective_weight_of_travelcost+region_pool[i][3]*alpha_coeff
        k=region_pool[i][4]
        #c2=sum(nodedik[x][k] for x in region_pool[i][0])
        #if abs(c2-cost)>0.001: print "check...", i,k,c2,cost,region_pool[i][1],region_pool[i][3]
        if fixed_cost_obj==1: ##[ulist,dis,demand,overload,k]
            cost+=facilityCost[k]
        if envy_service_objective==1:
            for j in region_pool[i][0]:
                if nodedij[j][k]>envy_service_distance:
                    d=nodedij[j][k]-envy_service_distance
                    cost+= nodes[j][3]* envy_objective_weight*d*d
        if envy_service_objective==4:
            for j in region_pool[i][0]:
                if nodedij[j][k]>envy_service_distance:
                    d=nodedij[j][k]-envy_service_distance
                    cost+= nodes[j][3]* envy_objective_weight*d
        costs.append(cost)
        #cost_stat.append(cost/region_pool[i][2])
        #cost_sum[k][0]+=1
        #cost_sum[k][1]+=cost/region_pool[i][2]
    #for k in range(num_districts):
    #    if cost_sum[k][0]>0:
    #        cost_sum[k][2]=cost_sum[k][1]/cost_sum[k][0]

    obj=""
    for i in range(len(region_pool)):
        obj+=costs[i]*variables[i]
    prob+=obj
    #for i in range(len(region_pool)):
    #    k=region_pool[i][4]
    #    if cost_stat[i]>cost_sum[k][2]:
    #        prob+= variables[i]==0
    #print region_pool
    #print costs
    rlist=[[] for i in range(num_units)]
    for j in range(len(region_pool)):
        for x in region_pool[j][0]:
            rlist[x].append(j)
    for i in range(num_units):
        s=""
        for x in rlist[i]:
            s+=variables[x]
        if spatial_contiguity==0:
            prob+=s >= 1
        else:
            prob+=s == 1
    if spatial_contiguity==0:
        for k in range(num_districts):
            s=""
            for i in range(len(region_pool)):
                if region_pool[i][4]!=k: continue
                s+=variables[i]
            if len(s)>0: prob+=s <= 1

    if adaptive_number_of_facilities==0:
        s=""
        for i in range(len(variables)):
            s+=variables[i]
        prob+= s==max_num_facility
    prob.writeLP("_spp.lp")
    #mip_mst_file=tempfile.mkstemp()[1].split("\\")[-1]

    vars=spp_mst()
    for x,v in vars: variables[x].setInitialValue(v)
    solver=0
    if mip_solver=='cbc': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.COIN_CMD(timeLimit=maxtime,mip=1,msg=solver_message,gapRel=mipgap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.cplex_api.CPLEX_CMD(msg=solver_message,warmStart=True,timeLimit=maxtime,options=['set parallel -1','set mip tolerances mipgap ' + str(mipgap)])
    if mip_solver=='gurobi': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.GUROBI_CMD(msg=solver_message,warmStart=True,options=[("MIPGap",mipgap),("TimeLimit",maxtime)])
    solver.setTmpDir() #=mip_file_path
    solver.actualSolve(prob)
    #if os.path.isfile(mip_mst_file): os.remove(mip_mst_file)
    if prob.status<=0:
       print "prob.status<0..."
       return prob.status #failer
    node_groups=[-1 for x in range(num_units)]
    centersID=[-1 for x in range(num_districts)]
    for v in prob.variables():
        if (v.varValue >= 0.99):
            items=v.name.split('_')
            i=int(items[1])
            k=region_pool[i][4]
            #print k,costs[i],facilityCost[k]
            centersID[k]=facilityCandidate[k]
            for x in region_pool[i][0]:
                node_groups[x]=k
    a=[ x for x in centersID if x>=0]
    print "spp locs:",len(a),a
    update_district_info()
    #for i in range(num_districts): 
    #    if district_info[i][1] >0: print i,district_info[i],facilityCost[i]
    #print 
    return 1 #success
        
def pop_selection(population):
    population1=copy.deepcopy(population)
    population1.sort(key=lambda x:x[0])
    population2=[] #delete identical solution
    #sol=[best_biobjective_global,best_centersID_global[:],best_solution_global[:],best_objective_global,best_objective_fcost_global,best_overload_global,0]
    #population2.append(copy.deepcopy(sol))
    population2.append(copy.deepcopy(population1[0]))
    sdiff=1
    if location_problem==3:
        sdiff=max_num_facility*5/100
        if sdiff<3: sdiff=3
        
    for x in population1:
        issimilar=0
        for y in population2:
            rlist=[i for i in range(num_districts) if x[1][i] != y[1][i]]
            if len(rlist)>=sdiff: continue
            else:
                if location_problem>=1:
                    issimilar=1
                    break
            ulist=[i for i in range(num_units) if x[2][i] != y[2][i]]
            #diffpop=sum(nodes[i][3] for i in ulist)
            #if len(ulist)<min(num_units*1.0/num_districts,num_units/30.0) and diffpop*100.0/total_pop < min(3.0,100.0/num_districts): #100.0/num_districts: #<5%
            #print len(ulist),
            if len(ulist)<num_units*(solution_similarity_limit/100.0):
                issimilar=1
                break
        if issimilar==0:
            population2.append(copy.deepcopy(x))
        #if len(population2)>=min(multi_start_count*3,10):
        #    break
    return population2
def cflpr_scp_model(maxr2,maxr,cover_pecentage,maxtime,mipgap): #bug, the covering constraint
    global node_groups
    global centersID
    if len(region_pool)<=10:
        print "len(region_pool)<=10", len(region_pool)
        return 0
    penalty=100000000
    prob = pulp.LpProblem("sdp_spp",pulp.LpMinimize)
    variables=[]
    costs=[]
    envyList=[]
    for i in range(len(region_pool)):
        x=pulp.LpVariable("x_" +str(i), 0, 1,pulp.LpBinary)
        variables.append(x)
        cost=region_pool[i][1]*envy_objective_weight_of_travelcost+region_pool[i][3]*penalty 
        k=region_pool[i][4]
        if fixed_cost_obj==1: 
            cost+=facilityCost[k]
        if envy_service_objective==1:
            envyobj=0.0
            for j in region_pool[i][0]:
                d=nodedij[j][k]
                if d>envy_service_distance:
                    envyobj+=nodes[j][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
            cost+=envyobj
        if envy_service_objective==4:
            envyobj=0.0
            for j in region_pool[i][0]:
                d=nodedij[j][k]
                if d>envy_service_distance:
                    envyobj+=nodes[j][3]*envy_objective_weight*(d-envy_service_distance)
            cost+=envyobj
        costs.append(cost)
        if envy_service_objective==2:
            envyobj=0.0
            for j in region_pool[i][0]:
                d=nodedij[j][k]
                if d>envy_service_distance:
                    envyobj+=nodes[j][3]*(d-envy_service_distance)*(d-envy_service_distance)
            envyList.append(envyobj)
    zvariables=[]
    for i in range(num_units):
        x=pulp.LpVariable("z_" +str(i), 0, 1,pulp.LpBinary)
        zvariables.append(x)
    obj=""
    for i in range(len(region_pool)):
        obj+=costs[i]*variables[i]
    prob+=obj

    for i in range(num_units):
        s=""
        for idx in range(len(region_pool)):
            if i in region_pool[idx][0]: 
                s+=variables[idx]
        if spatial_contiguity==0:
            prob+=s >= 1
        else:
            prob+=s == 1

    for i in range(num_units):
        s=""
        for idx in range(len(region_pool)):
            if i in region_pool[idx][0]: 
                k=region_pool[idx][4]
                if nodedij[i][k]<=maxr:s+=variables[idx]
        if spatial_contiguity==0:
            prob+=s- zvariables[i]>= 0
        else:
            prob+=s -zvariables[i] == 0

    for i in range(num_units):
        s=""
        for idx in range(len(region_pool)):
            if i in region_pool[idx][0]: 
                k=region_pool[idx][4]
                if nodedij[i][k]<=maxr2:
                    s+=variables[idx]
        prob+= s >= 1

    if spatial_contiguity==0:
        for k in range(num_districts):
            s=""
            for i in range(len(region_pool)):
                if region_pool[i][4]!=k: continue
                s+=variables[i]
            if len(s)>1: prob+=s <= 1

    for k in range(num_districts):
        if k not in facility_inclusion_list: continue
        s=""
        for i in range(len(region_pool)):
            if region_pool[i][4]!=k: continue
            s+=variables[i]
        if len(s)>0: prob+=s == 1
    s=""
    if adaptive_number_of_facilities==0:
        for i in range(len(region_pool)):
            s+=variables[i]
        prob+=s == max_num_facility
    s=""
    for i in range(num_units):
        s+=nodes[i][3]*zvariables[i]
    maxc=total_pop*cover_pecentage/100
    prob += s>=maxc

    if envy_service_objective==2:
        s=""
        for i in range(len(region_pool)):
            s += envyList[i]*variables[i]
        if len(s)>1:
            var=envy_coefficient_variation*envy_service_distance
            prob += s <= 0.5* total_pop* var*var
    prob.writeLP("_sclp.lp")
    #mip_mst_file=tempfile.mkstemp()[1].split("\\")[-1]

    vars=spp_mst()
    for x,v in vars: variables[x].setInitialValue(v)
    solver=0
    if mip_solver=='cbc': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.COIN_CMD(timeLimit=maxtime,mip=1,msg=solver_message,gapRel=mipgap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.cplex_api.CPLEX_CMD(msg=solver_message,warmStart=True,timeLimit=maxtime,options=['set parallel -1','set mip tolerances mipgap ' + str(mipgap)])
    if mip_solver=='gurobi': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.GUROBI_CMD(msg=solver_message,warmStart=True,options=[("MIPGap",mipgap),("TimeLimit",maxtime)])
    solver.setTmpDir() #=mip_file_path
    solver.actualSolve(prob)
    #if os.path.isfile(mip_mst_file): os.remove(mip_mst_file)
    if prob.status<=0:
       print "no solution! prob.status<=0..."
       return prob.status #failer
    node_groups=[-1 for x in range(num_units)]
    centersID=[-1 for x in range(num_districts)]
    for v in prob.variables():
        if (v.varValue >= 0.99):
            items=v.name.split('_')
            if items[0]!="x": continue
            i=int(items[1])
            k=region_pool[i][4]
            #print k,costs[i],facilityCost[k]
            centersID[k]=facilityCandidate[k]
            for x in region_pool[i][0]:
                node_groups[x]=k
    #a=[ x for x in centersID if x>=0]
    #print "spp locs:",len(a),a
    update_district_info()
    #for i in range(num_districts): 
    #    if district_info[i][1] >0: print i,district_info[i],facilityCost[i]
    #print 
    return 1 #success
def update_centers():
    global node_groups
    global centersID
    global time_update_centers
    if location_problem==1 or location_problem==0: return
    t=time.time()
    obj=biobjective
    sol=[-1 for x in range(num_units)]
    centers=[]
    for k in range(num_districts):
        if centersID[k]==-1: continue
        kn,ulist=update_center(k)
        for x in ulist: sol[x]=kn
        centers.append(kn)
        centersID[k]=-1
        #print [k,kn,k in ulist],
    node_groups=sol[:]
    for k in centers:
        centersID[k]=facilityCapacity[k]
    if location_problem==3:
        for i in range(num_units):
            for k in NearFacilityList[i]:
                if centersID[k]>=0:
                    node_groups[i]=k
                    break
    obj=biobjective
    update_district_info()
    update_best_solution()
    #print "updatecenters",biobjective-obj
    time_update_centers+=time.time()-t
    #print obj,obj-biobjective
def update_center(k):
    if location_problem==3: #pmp
        return update_center_pmp(k)
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    if ulist==[]: return k,[]
    best_cost=sum(nodedik[x][k] for x in ulist)
    best_center=k
    for i in range(num_districts):
            cost=sum(nodedik[x][i] for x in ulist)
            if cost<best_cost:
                best_cost=cost
                best_center=i
    return best_center,ulist
def update_center_pmp(k): #need debug for PMP with few cands
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    if ulist==[]: return k,[]
    best_cost=MAXNUMBER
    best_center=-1
    avgd=biobjective/num_units
    for i in range(num_districts):
        #if nodedik[i][k]>avgd/2: continue
        cost=sum(nodedik[x][i] for x in ulist)
        if cost<best_cost:
            best_cost=cost
            best_center=i
    return best_center,ulist

def location_check(key):
    if -1 in node_groups:
        arcpy_print("debug: "+str(key)+ " unit(s) not assigned! ")
        #return -1
    rlist=list(set(node_groups))
    if -1 in rlist: 
        arcpy_print("debug: "+str(key)+ " demand not assigned"+str(rlist))
        arcpy_print(str(node_groups))
        rlist.remove(-1)
    if len(rlist)>max_num_facility and adaptive_number_of_facilities==0:
        arcpy_print("debug: "+str(key)+ " too many facilities"+str(rlist))
        #return -1
    for k in range(num_districts):
        if k in rlist and centersID[k]==-1:
            arcpy_print("debug: "+str(key)+ " facilitiy not selected but used")
            #return -1
        if centersID[k]>=0 and k not in rlist:
            arcpy_print("debug: "+str(key)+ " facilitiy selected but not used")
            print k, district_info[k]
            print [x for x in centersID if x>=0]

            #return -1
        uid=centersID[k]
        if spatial_contiguity==1 and uid>-1 and node_groups[uid]!=k:
            arcpy_print("debug: "+str(key)+ " facilitiy unit assigned to other facility"+str(k))
            print k,uid, node_groups[uid]
            #return -1
    #return 1
   
def print_solution():
    arcpy_print("_______________final solution_______________")
    for i in range(num_units):
        s=""
        for x in nodes[i]:
            s+=str(x)+"\t"
        k=node_groups[i]
        kunit=centersID[k]
        s+=str(nodes[kunit][4])+"\t"
        selected=-1
        if i in facilityCandidate:
            selected=0
            if i in centersID:
                selected=1
        s+=str(selected)
        arcpy_print(s)
# return k centers by k-means
def k_means():
    global centersID
    global node_groups
    centers=[]
    centersxy=[[0.0,0.0] for x in range(max_num_facility)]
    sol=[-1 for x in range(num_units)]
    #random centers
    #for i in range(num_units):
    #    if nodes[i][3]>=avg_pop:
    #        centers.append(i)
    while 1:
        nid=random.randint(0,num_units-1)
        if nid not in centers:
            centers.append(nid)
        if len(centers)==max_num_facility:
            break
    for k in range(max_num_facility):
        cid=centers[k]
        centersxy[k][0]=nodes[cid][1]
        centersxy[k][1]=nodes[cid][2]
    loop=0
    distance_obj=MAXNUMBER
    #k-means
    while 1:
        loop+=1
        #print  centersxy
        sol = [-1 for x in range(num_units)]
        #random.shuffle(nodelist)
        total_d=0.0
        #assign node to center
        for i in range(num_units):
            cid=-1
            d=MAXNUMBER
            for k in range(max_num_facility):
                dx=nodes[i][1]-centersxy[k][0]
                dy=nodes[i][2]-centersxy[k][1]
                #dxy=pow(dx*dx+dy*dy,0.5)
                dxy=dx*dx+dy*dy
                if dxy<d:
                    d=dxy
                    cid=k
            sol[i]=cid
        for k in range(max_num_facility):
            klist=[x for x in range(num_units) if sol[x]==k]
            if len(klist)==0:
                continue
            dx=sum(nodes[x][1]*nodes[x][3] for x in klist)
            dy=sum(nodes[x][2]*nodes[x][3] for x in klist)
            dsum=sum(nodes[x][3] for x in klist)
            centersxy[k][0]=dx/dsum
            centersxy[k][1]=dy/dsum
        obj=0.0
        for i in range(num_units):
            k=sol[i]
            dx=nodes[i][1]-centersxy[k][0]
            dy=nodes[i][2]-centersxy[k][1]
            dxy=pow(dx*dx+dy*dy,0.5) * nodes[i][3]
            obj+=dxy
        if obj<distance_obj:
            distance_obj=obj
        else:
            break
        #print loop, distance_obj,obj
    centers=[]
    for k in range(max_num_facility):
        ulist=[x for x in range(num_units) if sol[x]==k]
        kdis=MAXNUMBER
        kcenter=-1
        for x in ulist:
            tmp=0.0
            for y in ulist:
                #tmp+=nodedij[y][x]*nodes[y][3]*(10+random.random())/10
                tmp+=nodedik[y][x]*(9.5+random.random())/10
            if tmp<kdis:
                kdis=tmp
                kcenter=x
        centers.append(kcenter)
    #print "k-means",distance_obj,centers
    centersID=[-1 for x in range(num_districts)]
    for x in centers: centersID[x]=x
    node_groups=[centers[x] for x in sol]
    update_district_info()
    return centers
def random_centers():
    global centersID
    centers=[]
    while 1:
        if len(centers)==max_num_facility: break
        k=random.randint(0,num_districts-1)
        if k not in centers: centers.append(k)
    for x in range(num_districts):centersID[x]=-1
    for x in centers: centersID[x]=facilityCandidate[x]
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break
    update_district_info()
    return 1
        
# return k centers by k-means
def k_medoids():
    global centersID
    global node_groups
    centers=[]
    sol=[-1 for x in range(num_units)]
    sol2=[-1 for x in range(num_units)]
   
    #while 1:
    #    nid=random.randint(0,num_units-1)
    #    if nid not in centers:
    #        centers.append(nid)
    #    if len(centers)==max_num_facility:
    #        break
    random_centers()
    centers=[x for x in range(num_districts) if centersID[x] >=0]
    loop=0
    distance_obj=MAXNUMBER
    #k-means
    while 1:
        loop+=1
        for i in range(num_units):
            for k in NearFacilityList[i]:
                if k in centers:
                    sol[i]=k
                    break
        obj=0.0
        for k in range(max_num_facility):
            ulist=[x for x in range(num_units) if sol[x]==centers[k]]
            if len(ulist)==0: continue
            clist=range(num_districts)
            cid=-1
            mindis=MAXNUMBER
            for i in clist:
                dsum=sum(nodedik[x][i] for x in ulist)
                if dsum<mindis:
                    mindis=dsum
                    cid=i
            centers[k]=cid
            obj+=mindis
            for i in ulist: sol2[i]=cid
        sol=sol2[:]
        if obj<distance_obj:
            distance_obj=obj
        else:
            break
        #print loop, distance_obj,obj,centers
    #print "k-means",distance_obj,centers
    centersID=[-1 for x in range(num_districts)]
    for x in centers: centersID[x]=facilityCandidate[x]
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break
    update_district_info()
    #update_centers()
    return 1

def k_medoids_sampling(allulist,cids):
    centers=cids
    num=len(allulist)
    if centers==[]:
        while 1:
            k=random.randint(0,num_districts-1)
            if k not in centers:
                centers.append(k)
            if len(centers)==max_num_facility:
                break
    sol=[-1 for x in allulist]
    loop=0
    distance_obj=MAXNUMBER
    #k-means
    while 1:
        loop+=1
        for idx in range(len(allulist)):
            i=allulist[idx]
            for k in NearFacilityList[i]:
                if k in centers:
                    sol[idx]=k
                    break
        obj=0.0
        for k in centers:
            ulist=[allulist[x] for x in range(num) if sol[x]==k]
            if len(ulist)==0: continue
            mindis=sum(nodedik[x][k] for x in ulist)
            cid=k
            for i in range(num_districts):
                dsum=sum(nodedik[x][i] for x in ulist)
                if dsum<mindis:
                    mindis=dsum
                    cid=i
            idx=centers.index(k)
            centers[idx]=cid
            obj+=mindis
        #print loop,distance_obj,sol
        if obj<distance_obj:
            distance_obj=obj
        else:
            break
    centers=list(set(centers))
    return distance_obj,centers

#wrighted k-medoids
#start_num: number of starts
#loop_num: number of loops 
#connective: constraints on spatial contiguity? 1 for yes and 0 for no    

def assign_ruin_recreate(idx):
    global time_ruin_recreate
    #if location_problem==4: return -1,1
    t=time.time()
    if idx<0 and len(ruin_oprators)<1: return -1
    ruin_idx=idx
    if idx<0:
        r=random.randint(0,len(ruin_oprators)-1)
        ruin_idx=ruin_oprators[r]
    #ruin_idx=0~4, diversification
    if ruin_idx==9: #move a few locations for PMP, PDP
        r_r_perb_center_locations()
    time_ruin_recreate[ruin_idx]+=time.time()-t
    return ruin_idx,1
def r_r_perb_center_locations(): #for pmp only, debuging
    global node_groups
    global centersID
    old_centersID=centersID[:]
    while 1:
        u=random.randint(0,num_units-1)
        k=node_groups[u]
        if centersID[k]==-1: continue
        for nk in NearFacilityList[u]:
            if nk==k: continue
            if centersID[nk]<0:
                centersID[k]=-1
                centersID[nk]=facilityCandidate[nk]
                break
        diff=sum(1 for x in range (num_districts) if centersID[x]!=old_centersID[x])
        if diff>=2: break
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break
    update_district_info()

def check_pmp_sol(pos):
    for i in range(num_units):
        k=node_groups[i]
        if centersID[k]<0: print i,k,"centersID[k]<=0",pos
        for j in NearFacilityList[i]:
            if centersID[j]<0: continue
            if k!=j and nodedik[i][k]!=nodedik[i][j] and spatial_contiguity==0: 
                print i,k,nodedik[i][k],j,nodedik[i][j],"k is not the nearest facility",pos
            break
    print pos,"best=",best_biobjective_global,"obj=",biobjective,"conn=",check_solution_continuality_feasibility(node_groups),
    print check_current_solution_continuality_feasibility()
def ils_pmp(numf,multistarts,timelimit,myseed):
    global multi_start_count
    global seed
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global objective
    global biobjective
    global objective_overload
    global node_groups
    global centersID
    global district_info
    #global node_neighbors
    global move_count
    global region_pool
    global pool_index
    global is_spp_modeling
    global spp_loops
    global pop_dis_coeff
    global all_solutions
    global assignment_operators_selected
    global max_num_facility
    global heuristic_time_limit
    global large_facility_cost
    global initial_solution_method
    global avg_dis_min
    global spatial_contiguity
    global location_tabu_list
    global max_exclusion_list
    max_num_facility=numf
    if location_problem==0:
        max_num_facility=num_districts
    initialize_instance()
    heuristic_time_limit=timelimit
    all_solutions=[]
    multi_start_count=multistarts
    seed=myseed
    if seed<0:
        seed=random.randint(0,100)
    random.seed(seed)
    arcpy_print("seed: "+str(seed))
    region_pool=[]
    t=time.time()
    ga_time=0.0
    best_biobjective_global = MAXNUMBER
    best_biobjective = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[]
    if is_spp_modeling>=1:
        pool_index=[[] for x in range(num_districts)]
    t_ga=time.time()
    multi_start=0

    population=generate_initial_solution()
    if len(max_exclusion_list)<1:
        max_exclusion_list=[0.0 for x in range(num_districts)]
    not_improved_global=0
    improved_time=time.time()
    loop=-1
    while 1:
        loop+=1
        r=random.random()
        sidx = int(min(multi_start_count,len(population))* pow(r,2)*0.999)
        r=random.random()
        node_groups=population[sidx][2][:]
        centersID=population[sidx][1][:]
        update_district_info()
        current=" current: " +str(int(population[sidx][0]))+" -> "
        not_improved_global+=1
        old_ids=centersID[:]
        best_obj=best_biobjective_global
        cur_obj=biobjective
        op=0
        sta=0
        #check_pmp_sol(0)
        if random.random()>0.5:
            idx,sta=r_r_new_location(-1)
        else:
            assign_ruin_recreate(-1)
            op=1
        if spatial_contiguity>=1: VND_local_search()
        #update_centers()
        #    check_pmp_sol(4)
        #check_pmp_sol(1)
        update_district_info()
        update_best_solution()
        update_region_pool_all()
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,max(0,total_pop-objective_supply)])
        s=""
        if biobjective<best_obj:  #0.1%
            s="*"
            impp=((biobjective-best_obj)/biobjective)*1000 #-0.1%
            not_improved_global+=int(max_loops_solution_not_improved*impp)
            if not_improved_global<0: not_improved_global=0
        elif biobjective<cur_obj-0.000001: 
            s="#"
        else: s="-"
        s=str(op)+" "+s
        bnf=sum(1 for x in best_centersID_global if x>=0)
        s+="Search loop "+str(loop) + " best: " +str(bnf)+" "+str(int(best_biobjective_global))+" "+str(int(best_objective_global))+" "+str(int(best_overload_global))
        nf=sum(1 for x in centersID if x>=0)
        s+=current + str(nf)+" "+str(int(biobjective))+" "+str(int(objective))+" "+str(int(objective_overload))
        f=[x for x in range(num_districts) if max_exclusion_list[x] <=best_biobjective_global+0.00001]
        n=sum(1 for x in range(num_districts) if best_centersID_global[x]>=0 and x in f) 
        s+=" "+ str(n)
        surplus=int(objective_supply-total_pop)
        psize=int(total_pop/nf)
        s+= " Info: " +str(len(population))+" " +str(len(potential_facilities)) +" " + str(surplus)+"/"+str(psize)+" " +str(not_improved_global)+ " "+str(int(time.time()-t_ga))
        arcpy_print(s)
        if not_improved_global >= max_loops_solution_not_improved: break
        population.sort(key=lambda x:x[0])
        population=pop_selection(population)
        #while len(population)>multi_start_count*5:
        #    population.pop()
        all_solutions=population
        if time.time()-t_ga > heuristic_time_limit:  break
    #post procedure

    population.sort(key=lambda x:x[0])
    node_groups=best_solution_global[:] #population[0][2][:]
    centersID=best_centersID_global[:]
    update_district_info()
    ga_time=time.time()-t_ga
    print "Heuristic solution:",biobjective,objective_fcost,objective,objective_overload,ga_time
    t_spp=time.time()
    if is_spp_modeling>=1:
        arcpy_print("SPP modelling..."+str(len(region_pool)) )
        sppmodel(heuristic_time_limit,0.000001)
        update_district_info()
        update_centers()
        update_best_solution()
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        print "spp solution:",biobjective,objective_fcost,objective,objective_overload,time.time()-t_spp, time.time()-t_spp
    print "final solution:",biobjective,objective_fcost,objective,objective_overload
    population.sort(key=lambda x:x[0])
    all_solutions=population
    sta=check_solution_continuality_feasibility(node_groups)
    print "Areal continuality(1 yes, 0 no):",sta
    if spatial_contiguity==1:
        if sta==0:    return "infeasible final solution: continuality"
    ##noprint "time",time.time()-t
    search_stat()
    print "location_tabu_list",len(location_tabu_list)
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

def ils_pmp_envy(numf,multistarts,timelimit,myseed):
    global multi_start_count
    global seed
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global objective
    global biobjective
    global objective_overload
    global node_groups
    global centersID
    global district_info
    #global node_neighbors
    global move_count
    global region_pool
    global pool_index
    global is_spp_modeling
    global spp_loops
    global pop_dis_coeff
    global all_solutions
    global assignment_operators_selected
    global max_num_facility
    global heuristic_time_limit
    global large_facility_cost
    global initial_solution_method
    global avg_dis_min
    global spatial_contiguity
    global location_tabu_list
    global max_exclusion_list
    max_num_facility=numf
    initialize_instance()
    heuristic_time_limit=timelimit
    all_solutions=[]
    multi_start_count=multistarts
    seed=myseed
    if seed<0:
        seed=random.randint(0,100)
    random.seed(seed)
    arcpy_print("seed: "+str(seed))
    region_pool=[]
    t=time.time()
    ga_time=0.0
    best_biobjective_global = MAXNUMBER
    best_biobjective = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[]
    if is_spp_modeling>=1:
        pool_index=[[] for x in range(num_districts)]
    t_ga=time.time()
    for idx in range(multistarts):
        best_biobjective = MAXNUMBER
        ulist=[]
        sample_size=min(num_units,30*max_num_facility)
        while len(ulist)<sample_size: 
            u=random.randint(0,num_units-1)
            if u not in ulist: ulist.append(u)
        random.shuffle(ulist)
        obj,centers=k_medoids_sampling(ulist,[])
        if len(centers)!=numf: 
            print "error in generating initial solution..."
            continue
        centersID=[-1 for x in range(num_units)]
        
        for x in centers: centersID[x]=facilityCandidate[x]
        for i in range(num_units):
            for k in NearFacilityList[i]:
                if centersID[k]<0: continue
                node_groups[i]=k
                break
        nf2=sum(1 for x in centersID if x>=0)
        update_envy_district_info()
        update_best_solution()
        all_solutions.append([biobjective,centersID[:],node_groups[:]])
        nf=sum(1 for x in centersID if x>=0)
        print "init solution", idx, nf2,nf,biobjective,objective
    not_improved_global=0
    population=all_solutions
    improved_time=time.time()
    loop=-1
    while 1:
        loop+=1
        r=random.random()
        sidx = int(min(multi_start_count,len(population))* pow(r,1)*0.999)
        r=random.random()
        node_groups=population[sidx][2][:]
        centersID=population[sidx][1][:]
        update_envy_district_info()
        current=" current: " +str(int(population[sidx][0]))+" -> "
        not_improved_global+=1
        old_ids=centersID[:]
        best_obj=best_biobjective_global
        cur_obj=biobjective
        strength=int(numf*0.05+0.5)
        if strength<2: strength=2
        else: strength=random.randint(2,strength)
        for i in range(strength):
            assign_ruin_recreate(9)
        '''for i in range(imploop):
            sta=pmp_Whitaker()
            if sta!=1:break
        update_envy_district_info()
        update_region_pool_all()
        update_best_solution()'''

        imploop=6
        for i in range(imploop):
            sta=pmp_Whitaker()
            update_envy_district_info()
            update_best_solution()
            update_region_pool_all()
            if sta!=1:break
        s=""
        if biobjective<best_obj:  #0.1%
            #del population[sidx]
            s="*"
            impp=((biobjective-best_obj)/biobjective)*1000 #-0.1%
            not_improved_global+=int(max_loops_solution_not_improved*impp)
            if not_improved_global<0: not_improved_global=0
            #population.append([biobjective,centersID[:],node_groups[:]])
        elif biobjective<cur_obj: 
            #del population[sidx]
            s="#"
            #population.append([biobjective,centersID[:],node_groups[:]])
        else: 
            s="-" 
            #del population[sidx]
        population.append([biobjective,centersID[:],node_groups[:]])
        population.sort(key=lambda x:x[0])
        population=pop_selection(population)
        if len(population)>max(10, multi_start_count*2):
            population.pop()
        bnf=sum(1 for x in best_centersID_global if x>=0)
        s+="Search loop "+str(loop) + " best: " +str(bnf)+" "+str(int(best_biobjective_global))+" "+str(int(best_objective_global))
        nf=sum(1 for x in centersID if x>=0)
        s+=current + str(nf)+" "+str(int(biobjective))+" "+str(int(objective))
        s+= " Info: " +str(len(population))+" " +str(not_improved_global)+ " "+str(int(time.time()-t_ga))
        arcpy_print(s)
        if not_improved_global >= max_loops_solution_not_improved: break
        #while len(population)>multi_start_count*5:
        #    population.pop()
        all_solutions=population
        if time.time()-t_ga > heuristic_time_limit:  break
    #post procedure

    population.sort(key=lambda x:x[0])
    node_groups=best_solution_global[:] #population[0][2][:]
    centersID=best_centersID_global[:]
    update_envy_district_info()
    ga_time=time.time()-t_ga
    print "Heuristic solution:",biobjective,objective_fcost,objective,objective_overload,ga_time
    t_spp=time.time()
    if is_spp_modeling>=1:
        arcpy_print("SPP modelling..."+str(len(region_pool)) )
        sppmodel(heuristic_time_limit,0.000001)
        #update_district_info()
        update_envy_district_info()
        #update_centers()
        update_best_solution()
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        print "spp solution:",biobjective,objective_fcost,objective,objective_overload,time.time()-t_spp, time.time()-t_spp
    print "final solution:",biobjective,objective_fcost,objective,objective_overload
    population.sort(key=lambda x:x[0])
    all_solutions=population
    #sta=check_solution_continuality_feasibility(node_groups)
    #print "Areal continuality(1 yes, 0 no):",sta
    #if spatial_contiguity==1:
    #    if sta==0:    return "infeasible final solution: continuality"
    ##noprint "time",time.time()-t
    #search_stat()
    #print "location_tabu_list",len(location_tabu_list)
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]


def search_stat():
    arcpy_print("----------------search statistics----------------------")
    arcpy_print("one unit move, move and time: "+ str(count_op[0])+ ", " +str(time_op[0]) )
    arcpy_print("two unit move, move and time: "+ str(count_op[1])+ ", " +str(time_op[1]) )
    arcpy_print("three unit move, move and time: "+ str(count_op[2])+ ", " +str(time_op[2]) )
    arcpy_print("location swap time: "+ str(time_location[0]) )
    arcpy_print("location drop time: "+ str(time_location[1]) )
    arcpy_print("location add time: "+ str(time_location[2]) )
    arcpy_print("location add-drop time: "+ str(time_location[3]) )
    arcpy_print("location multi-exchange time: "+ str(time_location[4]) )
    arcpy_print("r_r_reselect_location_pmp time: "+ str(time_location[5]) )
    arcpy_print("location TB heur. time: "+ str(time_location[7]) )
    arcpy_print("TB_Whitaker time: "+str(time_Whitaker))
    arcpy_print("location PMP sub_mip time: "+ str(time_location[8]) )
    arcpy_print("location CFLP sub_mip time: "+ str(time_location[9]) )
    arcpy_print("location PMP TB time: "+ str(time_pmp_re_location) )
    arcpy_print("repair time: "+ str(time_repair) )
    arcpy_print("check edge unit time: "+str(time_check_edge_unit))
    arcpy_print("update_centers time: "+ str(time_update_centers) )
    arcpy_print("spp regions: "+ str(len(region_pool)) )
    arcpy_print("spp pooling time: "+ str(time_spp) )
    arcpy_print("connectivity check time: "+ str(time_check))
    arcpy_print("time for ruin and recreate: " + str(time_ruin_recreate))
    #if spatial_contiguity==1:
    #    sta=check_solution_continuality_feasibility(best_solution_global)
    #    arcpy_print("solution on continuality (0 no, 1 yes) : "+str(sta))
    arcpy_print("----------------end of search statistics----------------")

def sol_stat():
    sol={}
    sol["num_facility"]=sum(1 for x in centersID if x>=0)
    sol["objective"]=biobjective
    sol["TotalCost"]=objective+objective_fcost
    sol["DCost"]=objective
    sol["Fcost"]=objective_fcost
    sol["WeightedEnvyValue"]=objective_envy
    if objective_envy>0:
        sol["EnvyValue"]=objective_envy/envy_objective_weight
    else:
        envy=0.0
        for i in range(num_units):
            k=node_groups[i]
            d=nodedij[i][k]
            if d>envy_service_distance:
                envy+=nodes[i][3]*(d-envy_service_distance)*(d-envy_service_distance)
        sol["EnvyValue"]=envy
    sol["OverloadPenalty"]=objective_pentalty_on_overload
    sol["RMaxPenalty"]=objective_pentalty_on_rmax
    sol["RCoverPenalty"]=objective_pentalty_on_covering
    return sol
  


    
    


def print_equality_measures():
    #print "-----------equality measures---------------"
    maxd=0.0
    mind=MAXNUMBER
    maxdev=0.0
    avgd=objective/total_pop
    absdev=0.0
    stddev=0.0
    Theil=0.0
    schutz=0.0
    for i in range(num_units): 
        k=node_groups[i]
        dis=nodedij[i][k]
        w=nodes[i][3]
        absdev+=w*abs(dis-avgd)
        stddev+=w*(dis-avgd)*(dis-avgd)
        if dis>maxd: maxd=dis
        if dis<mind: mind=dis
        if abs(dis-avgd)>maxdev: maxdev=abs(dis-avgd)
        '''if dis>0:
            a=dis*math.log(dis)-avgd*math.log(avgd)
            Theil+=w*a*a
        else:
            a=avgd*math.log(avgd)
            Theil+=w*a*a'''
        if dis>0:
            a=math.log(dis/avgd)
            Theil+=w*dis*a
        schutz+=w*abs(dis-avgd)
    equality_measures={}
    equality_measures["Mean"]=avgd
    #print "Centre", maxd
    equality_measures["Centre"]= maxd
    #print "Range",maxd-mind
    equality_measures["Range"]=maxd-mind
    #print "MaxDev", maxdev
    equality_measures["MaxDev"]=maxdev
    #print "MeanDev",absdev/total_pop
    equality_measures["MeanDev"]=absdev/total_pop
    #print "StdDev",math.sqrt(stddev/total_pop) #d.num_units
    equality_measures["Varance"]=stddev/total_pop
    equality_measures["StdDev"]=math.sqrt(stddev/total_pop)
    #print "Theil", Theil/avgd/total_pop #d.num_units
    equality_measures["VC"]=equality_measures["StdDev"] /avgd
   
    equality_measures["Theil"]= Theil/avgd/total_pop
    Gini=0.0
    for i in range(num_units): 
        k=node_groups[i]
        d1=nodedij[i][k]
        w1=nodes[i][3]
        for j in range(num_units): 
           k=node_groups[j]
           d2=nodedij[j][k]
           w2=nodes[j][3]
           Gini+=w1*w2*abs(d1-d2)
    #print "Gini",Gini/total_pop/total_pop/2/avgd #Gini/d.num_units/d.num_units/2/avgd
    equality_measures["Gini"]=Gini/total_pop/total_pop/2/avgd
    equality_measures["Schutz"]=schutz/total_pop/2/avgd
    return equality_measures