# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 00:08:25 2023

@author: YGebretsadik
"""
import os
import shapefile
import numpy as np
from datetime import datetime, timedelta
#import json
import pandas as pd
import networkx as nx
#import matplotlib.pyplot as plt
from pyvis.network import Network
from scipy.optimize import linprog
from scipy.io import loadmat
#import math
#from scipy.sparse import csr_matrix
from scipy.io import savemat
import time


def tic():
    return time.perf_counter()

def toc(start_time):
    elapsed_time = time.perf_counter() - start_time
    return elapsed_time

def read_shapefile(sf_shape):
    fields = [x[0] for x in sf_shape.fields][1:]
    records = [y[:] for y in sf_shape.records()]
    #records = sf_shape.records()
    shps = [s.points for s in sf_shape.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df

def forlder_structure(base):
    folder = {}
    folder['base']   = base
    folder['GIS']    = folder['base'] + '\\Dataset\\GIS'
    folder['irrdem'] = folder['base'] + '\\Dataset\\Irrigation_demands'
    folder['mandI']  = folder['base'] + '\\Dataset\\M_I'
    folder['infras'] = folder['base'] + '\\Dataset\\infrastucture'
    folder['runoff'] = folder['base'] + '\\Dataset\\Runoff'
    folder['result'] = folder['base'] + '\\result'
    return folder
def Hydropower_water_demand(t):
    #calculate water requirement for Hydropower 
    yval =[]
    fx=yval[t]
    return 
       
def Hydropower_ztransfer_function(t):
    #calculate slop of the gamma function  
    yval =[]
    fx=yval[t]
    return 

def operation(t):
    #for resOperations in OperationsClass
    #Hydropower_water_demand(t)
    yval=[0.96925615,0.968573625,0.969377357,0.9718157,1.00198023,0.98970285,1.078517061,1.113421202,1.004693486,0.987178867,0.981896792,0.975003639]
    yval=[1, 1, 1, 1, 1, 1, 1, 1 ,1 ,1 ,1,1,1,1,1,1,1]
    fx=yval[t]
    return fx 

def weird_division(n, d):
    return [(n/dd if dd else 0)  for  dd in d]  

def find_indices(larger_set, subset):
    indices = []
    for element in subset:
        try:
            index = larger_set.index(element)
            indices.append(index)
        except ValueError:
            pass
    return indices

def build_topology():
    #load shapefile into a dataframe and builds the topology
    return null 
    
if __name__ == "__main__":
    #%%  READ INPUT FILES into a dataframe  
    print("\033[H\033[J") 
    print('Loading files..')
    folder = forlder_structure('C:\\Users\\YGebretsadik\\OneDrive - CGIAR\\Documents\\experiment\\WorldBank')
    pNodedf = read_shapefile(shapefile.Reader(os.path.join(folder['GIS'], 'M_003_Nodes.shp')))
    pLinkdf = read_shapefile(shapefile.Reader(os.path.join(folder['GIS'], 'M_003_Links.shp')))
    Nu_nodes = pNodedf.shape[0]
    Nu_Links = pLinkdf.shape[0]
    Nu_Vars = Nu_Links+Nu_nodes
      
    #Remapping of Objects and Namse 
    pObjectNameMap = {'link' : 'link', 'res' : 'res', 'out' : 'out', 'cat' : 'cat', 'out_flow' : 'out_flow', 'end_Node' : 'source'}
    pNodedf['type'] = pNodedf['type'].map(pObjectNameMap)
   #%%  Create templet for display 
    NodeLooks = pd.DataFrame({'name' :['cat', 'link'    , 'out' , 'out_flow', 'res','source'],
                              'color':['#00FF00'  , '#808080' , '#FFA500',  '#000000', '#0000FF' ,'#FFFF00' ],
                              'size' :[30         , 1        , 20       ,          20,        30 ,30 ],
                              'shape':['circle'   ,'circle'   ,'box'     ,'star'      ,'triangle','hexagon']}) #https://graphviz.org/doc/info/shapes.html
    #%%  Get the Topology
    p1=np.array([[pLinkdf.iloc[x]['coords'][0][0] for x in range(Nu_Links)],[pLinkdf.iloc[x]['coords'][0][1] for x in range(Nu_Links)]]).T
    p2=np.array([[pLinkdf.iloc[x]['coords'][-1][0] for x in range(Nu_Links)],[pLinkdf.iloc[x]['coords'][-1][1] for x in range(Nu_Links)]]).T
    nl=np.array([[pNodedf.iloc[x]['coords'][0][0] for x in range(Nu_nodes)],[pNodedf.iloc[x]['coords'][0][1] for x in range(Nu_nodes)]]).T
    pfrom=np.array(pNodedf['ObjName'][np.array([np.argmin(((p1[x,0]-nl[:,0]) **2 + (p1[x,1]-nl[:,1])**2)**(0.5)) for x in range(Nu_Links)])])
    pto=np.array(pNodedf['ObjName'][np.array([np.argmin(((p2[x,0]-nl[:,0]) **2 + (p2[x,1]-nl[:,1])**2)**(0.5)) for x in range(Nu_Links)])])
    pTopology = np.column_stack((pfrom, pto))
   #%% Formulate the LP problem ====================================
    out_demand = loadmat(os.path.join(folder['irrdem'], 'Demands_low.mat'),simplify_cells=1)
    res_data=loadmat(os.path.join(folder['infras'], 'Dams_reser_con.mat'),simplify_cells=1)
    m_i_dem = loadmat(os.path.join(folder['mandI'], 'Demands.mat'),simplify_cells=1)
    cat_runoff=loadmat(os.path.join(folder['runoff'], 'Runoff_PF.mat'),simplify_cells=1)
    
    A1= np.int32(np.diag(np.ones(Nu_nodes)))
    A2= np.array([[(pNodedf['ObjName'][y] == pTopology[x-Nu_nodes,0])-(pNodedf['ObjName'][y] == pTopology[x-Nu_nodes,1])  for y in range(Nu_nodes)] for x in range(Nu_nodes,Nu_Links+Nu_nodes)])
    A=np.vstack(((np.hstack((A1,A2.T))),np.hstack((A1,np.int32(np.zeros(A2.shape)).T))))
    
    #A = (np.hstack((A1,A2.T)))
    indy_irr=[(pNodedf['ObjName'][y][0:6] == 'Irr_C_') for y in range(Nu_nodes)]
    indy_mni=[(pNodedf['ObjName'][y][0:4] == 'MNI_') for y in range(Nu_nodes)]
    indy_cat=[(pNodedf['type'][y] == 'cat') for y in range(Nu_nodes)]
    indy_res=[(pNodedf['type'][y] == 'res') for y in range(Nu_nodes)]
    indy_free_end=[(pNodedf['type'][y] == 'out_flow') for y in range(Nu_nodes)]
    indy_lnk=[(pNodedf['type'][y] == 'link') for y in range(Nu_nodes)]
    
    
    indx_cat=find_indices(cat_runoff['Cat_Names'].tolist(),pNodedf['ObjName'][indy_cat].tolist())
    indx_dem=find_indices(out_demand['pDemand'][:,0].tolist(),pNodedf['ObjName'][indy_irr].tolist())
    indx_res=find_indices(res_data['resNames'].tolist(),pNodedf['ObjName'][indy_res].tolist())
    indx_mni=find_indices(m_i_dem['pDemand'][:,0].tolist(),pNodedf['ObjName'][indy_mni].tolist())
    indx_hpflow=find_indices(pTopology[:,0].tolist(),pNodedf['ObjName'][indy_res].tolist())
    
    agents =np.array(pNodedf['ObjName'].tolist()+[a[0]+'-'+a[1] for a in pTopology])
    #%% Define scenrios 
    PrioritySCnearios={'Irrigation':[indy_mni,indy_irr,indy_res,indy_free_end],'Hydro':[indy_mni,indy_res,indy_irr,indy_free_end]}
    #RunoffSCnearios=[0,1,2,3,4,5,14]
    RunoffSCnearios=[0,1,2,3,4,5,14]
    StorageSCnearios={'Con':loadmat(os.path.join(folder['infras'], 'Dams_reser_con.mat'),simplify_cells=1),
                      'Reform':loadmat(os.path.join(folder['infras'], 'Dams_reser_reform.mat'),simplify_cells=1)}
    IrrigationSCenrios={'low':loadmat(os.path.join(folder['irrdem'], 'Demands_low.mat'),simplify_cells=1),
                       'High':loadmat(os.path.join(folder['irrdem'], 'Demands_high.mat'),simplify_cells=1),
                       'Current':loadmat(os.path.join(folder['irrdem'], 'Demands_current.mat'),simplify_cells=1)}
    
    
    PrioritySCnearios={'Irrigation':[indy_mni,indy_irr,indy_res,indy_free_end],'Hydro':[indy_mni,indy_res,indy_irr,indy_free_end]}
    RunoffSCnearios   =[0,1,2,3,4,5,14]
    StorageSCnearios  ={'Reform':loadmat(os.path.join(folder['infras'], 'Dams_reser_reform.mat'),simplify_cells=1)}
    IrrigationSCenrios={'High':loadmat(os.path.join(folder['irrdem'], 'Demands_high.mat'),simplify_cells=1)}
    
    #%% Solve accross Months========================================
    psolver='highs'
    pOptions = {"disp": False,'dual_feasibility_tolerance':1e-6,'primal_feasibility_tolerance':1e-3,'simplex_dual_edge_weight_strategy':'dantzig'}
    for strgsc in StorageSCnearios:
        for irrisc in IrrigationSCenrios:
            for Prioritysc in PrioritySCnearios:
                for runoffsc in RunoffSCnearios:
                    #%% solve
                    start_time = tic()
                    Res= np.zeros((Nu_Vars,360))
                    print('Simulation Started:'+datetime.fromtimestamp(start_time).strftime('%Y-%m-%d %H:%M:%S')+' Priority:'+Prioritysc)
                    #pPrio=[indy_mni,indy_irr,indy_res,indy_free_end]
                    #pPrio=[indy_mni,indy_res,indy_irr,indy_free_end]
                    pPrio=PrioritySCnearios[Prioritysc]
                    for t in range(360):
                        #initialize 
                        c=np.zeros(Nu_Vars)
                        b=np.zeros(2*Nu_nodes)
                        xbounds=np.array([(0.0,0.0)]*Nu_nodes + [(0.0,99.9e9)]*Nu_Links)
                       
                        b[0:Nu_nodes][indy_cat] = abs(cat_runoff['Runoff'][indx_cat,runoffsc,t]) 
                        b[0:Nu_nodes][indy_res] = StorageSCnearios[strgsc]['Resdata'][indx_res,1]*1e6 if t==0 else Res[0:Nu_nodes,t-1][indy_res] 
                        pWanted=np.zeros(Nu_nodes);
                        
                        pWanted[indy_irr]=IrrigationSCenrios[irrisc]['pDemand'][indx_dem,t%12+1]*1e6;
                        pWanted[indy_res]=operation(t%12)*(StorageSCnearios[strgsc]['Resdata'][indx_res,1]*1e6 if t==0 else Res[0:Nu_nodes,t-1][indy_res])
                        pWanted[indy_mni]=m_i_dem['pDemand'][indx_mni,t%12+1]*1e6
                        pWanted[indy_free_end]=99.9e9
                        ind_hi_cer=np.arange(Nu_nodes)
                        ind_lo_cer=np.arange(Nu_nodes,2*Nu_nodes)
                        for prio in pPrio:
                            evaluation_index=prio
                            c=np.array([1]*Nu_nodes + [0]*Nu_Links)
                            c[0:Nu_nodes][evaluation_index] = -1
                            b[Nu_nodes:2*Nu_nodes][evaluation_index] = pWanted[evaluation_index]
                            OptRes=linprog(c, A_ub=A[ind_lo_cer[evaluation_index],:], b_ub=b[ind_lo_cer[evaluation_index]], A_eq=A[ind_hi_cer,:], b_eq=b[ind_hi_cer],bounds=(0.0,None),method= psolver,options = pOptions)
                            Res[:,t]= abs(OptRes['x'])
                            b[Nu_nodes:2*Nu_nodes][evaluation_index]=Res[0:Nu_nodes,t][evaluation_index]
                            ind_hi_cer=np.concatenate((ind_hi_cer,ind_lo_cer[evaluation_index])) 
                        print('time step->' + str(t+1)+',scenario:'+ str(runoffsc))
                    #%% Post Process
                    #FLowatHP=np.array( [(ResultMCM[(A[0:Nu_nodes,:][indy_res,:]==1)[x,:],:]) for x in range(len(indy_res))] ).squeeze()
                    ResultMCM=abs(Res/1e6)
                    nethead=StorageSCnearios[strgsc]['Resdata'][indx_res,7]
                    installedcap=StorageSCnearios[strgsc]['Resdata'][indx_res,5]
                    pHp_flowm3=(Res[Nu_nodes:2*Nu_nodes,:][indx_hpflow])/(3600)
                    
                    pEnergy=np.array([pHp_flowm3[:,t]*1000*nethead*0.85*0.95*0.93*9.806 for t in range(360)] )/1e9 #in GWh
                    HP_res ={'Names':agents[0:Nu_nodes][indy_res],'value':2*pEnergy,'Install Cap':installedcap}
                    Irr_res={'Names':agents[0:Nu_nodes][indy_irr],'value':2*ResultMCM[0:Nu_nodes,:][indy_irr,:]}
                    Flow_res={'Names':agents[Nu_nodes:2*Nu_nodes],'value':ResultMCM[Nu_nodes:2*Nu_nodes,:]}
                    Res_res={'Names':agents[0:Nu_nodes][indy_res],'value':ResultMCM[0:Nu_nodes,:][indy_res,:]}
                    MnI_res={'Names':agents[0:Nu_nodes][indy_mni],'value':ResultMCM[0:Nu_nodes,:][indy_mni,:]}
                    
                    a = np.arange(20)
                    mdic = {"a": a, "label": "experiment"}
                    SimulationResult={
                     'HydroPower':HP_res , 
                     'flow':Flow_res ,
                     'Reservoir':Res_res,
                     'Irrigation':Irr_res,  #calibration Unit 
                     'M_n_I':MnI_res  
                     }
                    savemat(os.path.join(folder['result'], 'Result_2'+str(strgsc)+'_'+str(irrisc)+'_'+ str(runoffsc)+'_'+Prioritysc +'.mat'), SimulationResult)
                    elapsed_time = toc(start_time)
                    print(f"Elapsed time: {elapsed_time} seconds")
    #%% Display network in pyvis ======================================
    plotindx=0
    pModel = nx.DiGraph()
    for x in range(Nu_nodes):
        pModel.add_node(pNodedf.iloc[x]['ObjName'], 
                        demand=5,
                        color=NodeLooks[NodeLooks['name'] ==pNodedf.iloc[x]['type']]['color'].iloc[0],
                        shape=NodeLooks[NodeLooks['name']==pNodedf.iloc[x]['type']]['shape'].iloc[0],
                        size= NodeLooks[NodeLooks['name']==pNodedf.iloc[x]['type']]['size'].iloc[0],
                        x=(pNodedf.iloc[x]['coords'][0][0]-40)*500,
                        y=500-(pNodedf.iloc[x]['coords'][0][1]+8)*500
                      )
    for x in range(Nu_Links):
        pModel.add_edge(pTopology[x,0],pTopology[x,1],
                        weight=ResultMCM[Nu_nodes+x,plotindx]/30, 
                        capacity=5,
                        color = '#0789f4',
                        label="(" + "{:.3f}".format(np.mean(ResultMCM[Nu_nodes+x,:])*12) + ")",
                        font={'color': 'red','bold': True, 'face': 'Cambria'},
                        arrows='to')
   
    nt = Network('1000px', '1900px',select_menu=True,)
    nt.from_nx(pModel)
    nt.prep_notebook()
    nt.show_buttons(filter_=True)
    nt.toggle_physics(False)
    #nt.force_atlas_2based(overlap= 1)
    nt.show('nx2.html')
     
    
   