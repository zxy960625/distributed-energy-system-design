# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

# Parameters
building_load_design=pd.read_excel('parameter_setting/demand_side.xlsx')
weather=pd.read_csv('parameter_setting/Hong_Kong_SAR_CH-hour.csv')
period=24
period1 = 24
model = gp.Model('DES_EM')
dt= pd.date_range("2015-01-01", periods=8760, freq="H")
building_load_design['date']=dt
weather['date']=dt
building_load_design=building_load_design.set_index('date')
building_load_design=building_load_design['2015-07-15']
weather=weather.set_index('date')
weather=weather['2015-07-15']
#design time interval
CL_design=building_load_design.iloc[:,0].values.ravel()
EL_design=building_load_design.iloc[:,1].values.ravel()
I=weather.iloc[:,0].values.ravel()
T=weather.iloc[:,1].values.ravel()
#price
price_guangdong=pd.read_excel('parameter_setting/E_Guangdong5.xlsx')
price_caiso=pd.read_excel('parameter_setting/energy_price.xlsx')
price_guangdong=price_guangdong[:24]
price_caiso=price_caiso[:24]
P_e=price_guangdong.iloc[:,1].values.ravel()
P_g=price_caiso['Gas'].values.ravel()
#gas turbin
#capacity
# cap_gt = model.addVar(vtype=GRB.CONTINUOUS, name="cap_gt")
cap_gt = 6000
#cost
c_gt=1100
#efficiency
eff_ele=0.42
eff_heat=0.38
#gastopower
gas_power=12.62
#amount of gas
gas_gt_power=model.addVars(period,vtype=GRB.CONTINUOUS, name="gas_gt_power")
gas_gt_heat=model.addVars(period,vtype=GRB.CONTINUOUS, name="gas_gt_heat")
gas_gt=model.addVars(period,vtype=GRB.CONTINUOUS, name="gas_gt")
#model power_gt_fuel
power_gt=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_gt")
#amount of electricity produced by DG

gt_power=model.addConstrs(power_gt[i] == eff_ele* gas_gt_power[i]*gas_power
                          for i in range(period))

#constraints
cons_gt = model.addConstrs(power_gt[i] <= cap_gt*2
                           for i in range(period))
#model heat
heat_gt=model.addVars(period,vtype=GRB.CONTINUOUS, name="heat_gt")
#amount of heat produced by DG
gt_heat=model.addConstrs(heat_gt[i] == eff_heat * gas_gt_heat[i]*gas_power
                          for i in range(period))
#capacity
# cap_ac = model.addVar(vtype=GRB.CONTINUOUS,name="cap_ac")
# cap_ec = model.addVar(vtype=GRB.CONTINUOUS, name="cap_ec")
max_cl=32560
design_cl=max_cl*1.1
# model.addConstr(cap_ac== cap_gt*0.38/0.42*1.2/3)
# model.addConstr(cap_ec== (design_cl-6*cap_ac)/6)
cap_ac= 2000
cap_ec= 4000
#cooling load
cl_ec=model.addVars(period,vtype=GRB.CONTINUOUS,name="cl_ec")
cl_ac=model.addVars(period,vtype=GRB.CONTINUOUS, name="cl_ac")
power_ec_single=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_ec_single")
cl_ec_single=model.addVars(period,vtype=GRB.CONTINUOUS, name="cl_ec_single")
r_ec=model.addVars(period,vtype=GRB.CONTINUOUS, name="r_ec")
n_ec=model.addVars(period,vtype=GRB.CONTINUOUS, name="n_ec")
n_act_ec=model.addVars(period,vtype=GRB.INTEGER, name="n_act_ec")
r=model.addVars(period,vtype=GRB.CONTINUOUS, name="r")
#capital cost
c_ec=151.5
c_ac=180
#ec&ac efficiency
cop_ac=1.2
g0=model.addVars(period,vtype=GRB.BINARY,name="g0")
g1=model.addVars(period,vtype=GRB.BINARY,name="g1")
g2=model.addVars(period,vtype=GRB.BINARY,name="g2")
M=10e6
#OPEN WHEN R_EC>0.2
for i in range(period):
    model.addConstr(g0[i]+g1[i]+g2[i]==1)
    model.addConstr(n_ec[i]==cl_ec[i]/cap_ec)
    model.addConstr(r[i]+n_act_ec[i]==n_ec[i]+1)
    model.addConstr(r[i]<=1)
    model.addConstr(cl_ec_single[i]*n_act_ec[i]==cl_ec[i])
    model.addConstr(cl_ec_single[i]<=0+M*(1-g0[i]))
    model.addConstr(-cl_ec_single[i] <= 0 + M*(1 - g0[i]))
    model.addConstr(power_ec_single[i]<=0+M*(1-g0[i]))
    model.addConstr(-power_ec_single[i] <= 0 + M * (1 - g0[i]))
    model.addConstr(cl_ec_single[i]<=0+M*(1-g1[i])+0.5*cap_ec)
    model.addConstr(-cl_ec_single[i] <= 0 + M * (1-g1[i])-0.2 * cap_ec)
    model.addConstr(power_ec_single[i]<=M*(1-g1[i])+cl_ec_single[i]/3)
    model.addConstr(-power_ec_single[i] <= 0 + M * (1 - g1[i])-cl_ec_single[i]/3)
    model.addConstr(cl_ec_single[i]<=0+M*(1-g2[i])+cap_ec)
    model.addConstr(-cl_ec_single[i] <= 0 + M * (1 - g2[i])-0.5*cap_ec)
    model.addConstr(power_ec_single[i]<=0+M*(1-g2[i])+cl_ec_single[i]/5)
    model.addConstr(-power_ec_single[i]<=0+M*(1 -g2[i])-cl_ec_single[i]/5)
    #model.addGenConstrMax(r_ec[i],[y[i],0.2],name="maxconstr")
    #model.addGenConstrPWL(r_ec[i],cop_ec[i],x_pw,ec_pw)
#constraints
cons_cl_ac = model.addConstrs(cl_ac[i]<= cap_ac*6
                         for i in range(period))
cons_cl_ec = model.addConstrs(cl_ec[i]<= cap_ec*6
                         for i in range(period))
#model absorption chiller load
ac_cl=model.addConstrs(cl_ac[i]== cop_ac * heat_gt[i]
                        for i in range(period))
#model electric chiller load
# ec_cl =model.addConstrs(cl_ec[i]== cop_actual_ec[i]* power_ec[i]
#                         for i in range(period))
#battery
#capacity
cap_bat = 10000
#cost
cha_max_bat=0.2
dis_max_bat=0.5
c_bat=179

#soc_bat
soc_bat=model.addVars(period,vtype=GRB.CONTINUOUS, name="soc_bat")
#power
power_bat_dis=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_bat_dis")
power_bat_cha=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_bat_cha")
#state
bat_cha= model.addVars(period,vtype=GRB.BINARY, name="bat_cha")
bat_dis = model.addVars(period,vtype=GRB.BINARY, name="bat_dis")
#constraints
model.addConstrs(bat_cha[i]+bat_dis[i]<=1
                         for i in range(period))

cons_bat_max = model.addConstrs(0.95*cap_bat>= soc_bat[i]
                         for i in range(period))
cons_bat_min = model.addConstrs(0.1*cap_bat<= soc_bat[i]
                         for i in range(period))
cons_bat_charge_max = model.addConstrs(cha_max_bat*cap_bat>= power_bat_cha[i]
                         for i in range(period))
cons_bat_discharge_max=model.addConstrs(dis_max_bat*cap_bat>= power_bat_dis[i]
                         for i in range(period))
#soc_bat
model.addConstr(soc_bat[0]==0.1*cap_bat)
model.addConstrs(soc_bat[i]==0.1*cap_bat
                 for i in range(period1-1,period,period1))
model.addConstrs((soc_bat[i+1]== soc_bat[i]+0.95*bat_cha[i]*power_bat_cha[i]-bat_dis[i]*power_bat_dis[i]/0.95
                  for i in range(period-1)))
model.addConstrs(soc_bat[i-1]>=power_bat_dis[i]
                 for i in range(1,period))
#cooling storage
cap_cs = 5000
#cost
cs_cha_max_bat=0.2
c_cs=51
#soc
soc_cs=model.addVars(period,vtype=GRB.CONTINUOUS, name="soc_cs")
#power
power_cs_dis=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_cs_dis")
power_cs_cha=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_cs_cha")
#state
cs_cha= model.addVars(period,vtype=GRB.BINARY, name="cs_cha")
cs_dis = model.addVars(period,vtype=GRB.BINARY, name="cs_dis")
#constraints
model.addConstrs(cs_cha[i]+cs_dis[i]<=1
                         for i in range(period))

cons_cs_max = model.addConstrs(0.9*cap_cs>= soc_cs[i]
                         for i in range(period))
cons_cs_max = model.addConstrs(0.3*cap_cs<= soc_cs[i]
                         for i in range(period))
cons_cs_charge_max = model.addConstrs(cs_cha_max_bat*cap_cs>= power_cs_cha[i]
                         for i in range(period))
cons_cs_discharge_max=model.addConstrs(cs_cha_max_bat*cap_cs>= power_cs_dis[i]
                         for i in range(period))
#soc
model.addConstr(soc_cs[0]==0.3*cap_cs)
model.addConstrs(soc_cs[i]==0.3*cap_cs
                 for i in range(period1-1,period,period1+24))
model.addConstrs((soc_cs[i+1]== 0.99*soc_cs[i]+cs_cha[i]*power_cs_cha[i]*0.88-cs_dis[i]*power_cs_dis[i]/0.88
                  for i in range(period-1)))
model.addConstrs(soc_cs[i-1]>=power_cs_dis[i]
                 for i in range(1,period))
#PV
#capacity
cap_pv = 157
# t_cell = model.addVars(period,vtype=GRB.CONTINUOUS, name="t_cell")
# model.addConstrs(t_cell[i]==T[i]+0.0256*I[i]
#                  for i in range(period))
#1 square meter= 157wp
#1Wp=$0.82
#cost
c_pv=180
#efficiency
eff_pv=0.15
#pv_power
# pv_power=model.addVars(period,vtype=GRB.CONTINUOUS, name="pv_power")
#amount of gas
pv_power=np.zeros(24)
t_cell=np.zeros(24)
for i in range(period):
    t_cell[i] == T[i] + 0.0256 * I[i]
    pv_power[i] = I[i]*eff_pv*cap_pv*(1-0.0037*t_cell[i]+0.0037*25)*0.001
#Grid
#power_gird
power_grid_buy=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_grid_buy")
power_grid_sell=model.addVars(period,vtype=GRB.CONTINUOUS, name="power_grid_sale")
#Demand response
dr_up=model.addVars(period,vtype=GRB.CONTINUOUS, name="dr_up")
dr_down=model.addVars(period,vtype=GRB.CONTINUOUS, name="dr_down")
dr_up_total=model.addVar(vtype=GRB.INTEGER, name="dr_up_total")
dr_down_total=model.addVar(vtype=GRB.INTEGER, name="dr_down_total")
dr_up_open=model.addVars(period,vtype=GRB.BINARY, name="dr_up_open")
dr_down_open=model.addVars(period,vtype=GRB.BINARY, name="dr_down_open")
model.addConstrs(dr_up_open[i]+dr_down_open[i]<=1 for i in range(period))
model.addConstrs(dr_up_open[i]*dr_up[i]<=0.1*EL_design[i] for i in range(period))
model.addConstrs(dr_down_open[i]*dr_down[i]<=0.1*EL_design[i] for i in range(period))
# model.addConstrs(gp.quicksum(dr_up_open[i]*dr_up[i] for i in range(k*period1,(k+1)*period1)) ==
#                    gp.quicksum(dr_down_open[i]*dr_down[i] for i in range(k*period1,(k+1)*period1)) for k in range(365))
model.addConstr(gp.quicksum(dr_up_open[i]*dr_up[i] for i in range(24)) ==gp.quicksum(dr_down_open[i]*dr_down[i] for i in range(24)))
#passive thermal storage
#fast regulation
sp_gt=model.addVars(period,vtype=GRB.CONTINUOUS, name="sp_gt")
rg_bat=model.addVars(period,vtype=GRB.CONTINUOUS, name="rg_down_bat")
#gt
model.addConstrs(power_gt[i]+sp_gt[i]<=cap_gt for i in range(period))
# model.addConstrs(rg_down_gt[i]==rg_up_gt[i] for i in range(period))
# #ec
# model.addConstrs(rg_ec[i]<=CL_design[i]/cop_ec*0.15*0.3 for i in range(period))#temperature control
# model.addConstrs(rg_ec[i]/0.3*+power_ec[i]<=cap_ec/cop_ec for i in range(period))
# model.addConstrs(rg_ec[i]<=power_grid_buy[i] for i in range(period))

#model.addConstrs(rg_up_ec[i]+sp_ec[i]<=0.2*power_ec[i] for i in range(period))
#model.addConstrs(rg_down_ec[i]+sp_ec[i]<=0.2*power_ec[i] for i in range(period))
# model.addConstrs(rg_ec[i]+power_ec[i]<=cap_ec for i in range(period))
#bat
model.addConstrs(soc_bat[i-1]+rg_bat[i]*0.95<= 0.95*cap_bat for i in range(1,period,1))
model.addConstrs(soc_bat[i-1]-rg_bat[i]/0.95>= 0.1*cap_bat for i in range(1,period,1))
model.addConstrs(rg_bat[i]<= cha_max_bat*cap_bat for i in range(period))
model.addConstrs(rg_bat[i]<=dis_max_bat*cap_bat for i in range(period))
model.addConstrs(soc_bat[i-1]>=rg_bat[i] for i in range(1,period,1))
#energy balance model
cl_distribution=model.addVars(period,vtype=GRB.CONTINUOUS, name="cl_distribution")
model.addConstrs(cl_distribution[i]==0.12*(CL_design[i]-cs_cha[i]*power_cs_cha[i]+cs_dis[i]*power_cs_dis[i])
                 for i in range(period))
energy_balance= model.addConstrs((bat_dis[i]*power_bat_dis[i]+ power_grid_buy[i]+ power_gt[i]+dr_down_open[i]*dr_down[i]+pv_power[i]
                                 ==bat_cha[i]*power_bat_cha[i]+n_act_ec[i]*power_ec_single[i]+EL_design[i]+dr_up_open[i]*dr_up[i]+power_grid_sell[i]+cl_distribution[i]
                                 for i in range(period)))

#cooling load model
cooling_balance = model.addConstrs(cl_ec[i]+cl_ac[i]+cs_dis[i]*power_cs_dis[i]-cs_cha[i]*power_cs_cha[i]==CL_design[i]
                                    for i in range(period))
#life cycle
r=0.07

year_gt=0.09439
year_ec=0.1098
year_ac=0.1098
year_bat=0.1424
year_pv=0.09439
year_cs=0.09439
e_gas=1.93
e_grid=0.804
c_carbon=0.05019*0.16
#objective functions

operation=gp.quicksum(P_g[i]*gas_gt[i]+P_e[i]*power_grid_buy[i]
                      for i in range(period))

carbon_emission=gp.quicksum(c_carbon*e_gas*gas_gt[i]+e_grid*c_carbon*power_grid_buy[i]
                            for i in range(period))
investment=year_gt*c_gt*2*cap_gt+ year_ec*c_ec*6*cap_ec+ year_ac*c_ac*6*cap_ac +year_bat*c_bat*cap_bat+year_pv*c_pv*cap_pv+year_pv*c_cs*cap_cs
P_sp=20*10e-4*0.16
P_rg=20*10e-4*0.16
as_revenue=gp.quicksum(P_sp*sp_gt[i]+
                        P_rg*rg_bat[i]
                        for i in range(period))
#maintenance cost
maintenance=0.03*investment
model.setObjective(investment+operation+maintenance+carbon_emission-as_revenue)
model.Params.TimeLimit = 10000
model.params.NonConvex = 2
model.optimize()
out_df_cap=pd.Series(['1','1','2','2','1','1','1'],
                     index=['object','cap_gt','cap_ac','cap_ec','cap_bat','cap_pv','cap_cs'])
out_df_ope=pd.DataFrame(columns=['cl_ec','cl_ac','CL_design','EL_design','power_cs_dis','power_cs_cha','soc_cs',
                                 'power_grid_buy','power_gt','dr_up','dr_down','power_bat_dis','power_ec','cl_distribution'
                                 'power_bat_cha','power_pv','soc_bat','gas_gt','e_price','rg_bat','sp_gt'])
                                 # 'rg_up_ec','rg_down_ec','rg_up_gt','rg_down_gt','rg_up_bat','rg_down_bat',

# out_df_cap['object']=model.ObjVal
# out_df_cap['cap_gt']=cap_gt.x
# out_df_cap['cap_ac']=cap_ac.x
# out_df_cap['cap_ec']=cap_ec.x
# out_df_cap['cap_bat']=cap_bat.x
# out_df_cap['cap_pv']=cap_pv.x
# out_df_cap['cap_cs']=cap_cs.x
out_df_ope['cl_distribution']=model.getAttr('X', cl_distribution.values())
out_df_ope['power_ec_single']=model.getAttr('X',power_ec_single.values())
out_df_ope['cl_ec']=model.getAttr('X', cl_ec.values())
out_df_ope['cl_ac']=model.getAttr('X', cl_ac.values())
out_df_ope['CL_design']=CL_design
out_df_ope['power_cs_dis']=model.getAttr('X', power_cs_dis.values())
out_df_ope['power_cs_cha']=model.getAttr('X', power_cs_cha.values())
out_df_ope['soc_cs']=model.getAttr('X',soc_cs.values())
out_df_ope['EL_design']=EL_design
out_df_ope['n_act_ec']=model.getAttr('X',n_act_ec.values())
out_df_ope['power_grid_buy']=model.getAttr('X',power_grid_buy.values())
out_df_ope['power_gt']=model.getAttr('X', power_gt.values())
out_df_ope['dr_up']=model.getAttr('X', dr_up.values())
out_df_ope['dr_down']=model.getAttr('X', dr_down.values())
out_df_ope['power_bat_dis']=model.getAttr('X', power_bat_dis.values())
out_df_ope['power_bat_cha']=model.getAttr('X',power_bat_cha.values())
out_df_ope['soc_bat']=model.getAttr('X',soc_bat.values())
# out_df_ope['power_pv']=model.getAttr('X',pv_power.values())
out_df_ope['gas_gt']=model.getAttr('X',gas_gt.values())
out_df_ope['sp_gt']=model.getAttr('X',sp_gt.values())
out_df_ope['rg_bat']=model.getAttr('X',rg_bat.values())
# out_df_ope['rg_up_gt']=model.getAttr('X',rg_up_gt.values())
# out_df_ope['rg_down_gt']=model.getAttr('X',rg_down_gt.values())
# out_df_ope['rg_up_ec']=model.getAttr('X',rg_up_ec.values())
# out_df_ope['rg_down_ec']=model.getAttr('X',rg_up_ec.values())
# out_df_ope['rg_up_bat']=model.getAttr('X',rg_up_ec.values())
# out_df_ope['rg_down_bat']=model.getAttr('X',rg_up_ec.values())
out_df_ope['e_price']=P_e
start=pd.Timestamp(2015,1,1,0)
end=pd.Timestamp(2015,1,1,23)
dt1=pd.date_range(start,end,freq='H')
out_df_ope=out_df_ope.set_index(dt1)
# out_df_cap.to_excel('D:/PycharmProjects/pythonProject/out/'+'e10-Cap-trail.xlsx')
out_df_ope.to_excel('D:/PycharmProjects/pythonProject/out/a'+'trail.xlsx')
