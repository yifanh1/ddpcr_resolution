#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from brokenaxes import brokenaxes
plt.style.use('seaborn-white')


# In[3]:


#AT-TAQMAN
data=pd.read_csv("RT_PCR_TaqMan_annealing_temperature_probes.csv")
x=np.array(data['AT'])
y16=np.array(data['dCt (HPV16)'])
y18=np.array(data['dCt (HPV18)'])
fig=plt.figure(figsize=(12,12))#create the figure
#subset1
plt.subplot(2,1,1)
from scipy.interpolate import interp1d
xnew = np.linspace(x.min(),x.max(),300) 
plt.plot(x,y16,'.',color='royalblue')
plt.plot(x,y18,'.',color='maroon')
inter16 = interp1d(x, y16 , kind='cubic')
inter18 = interp1d(x, y18, kind='cubic')

y16_smooth = inter16(xnew)
y18_smooth=inter18(xnew)
plt.plot(xnew,y16_smooth,label='HPV16',alpha=0.7)
plt.plot(xnew,y18_smooth,label='HPV18',color='firebrick',alpha=0.7)


#set the title and x,y labels
plt.title('Taqman Annealing Temperature')
plt.xlabel('Annealing Tempearature ($^\circ$C)')
plt.ylabel('$\Delta$Ct ')
plt.legend()
#plt.gcf()--current figure,plt.gca()currenr axis

#adjust the axis limits
#plt.xlim(-1,1)
plt.ylim(14,21)
#plt.axis('tight')


# In[4]:


#AT-EVAGREEN
data=pd.read_csv("RT_PCR_evagreen_annealing_temperature_primers.csv")
x=np.array(data['AT'])
y16=np.array(data['dCtHPV16'])
y18=np.array(data['dCtHPV18'])
fig=plt.figure(figsize=(12,12))#create the figure
#subset2
plt.subplot(2,1,2)
from scipy.interpolate import interp1d
xnew = np.linspace(x.min(),x.max(),300) 
plt.plot(x,y16,'.',color='royalblue')
plt.plot(x,y18,'.',color='maroon')
inter16 = interp1d(x, y16 , kind='cubic')
inter18 = interp1d(x, y18, kind='cubic')

y16_smooth = inter16(xnew)
y18_smooth=inter18(xnew)
plt.plot(xnew,y16_smooth,label='HPV16',alpha=0.7)
plt.plot(xnew,y18_smooth,label='HPV18',color='firebrick',alpha=0.7)


#set the title and x,y labels
plt.title('Evagreen Annealing Temperature')
plt.xlabel('Annealing Tempearature ($^\circ$C)')
plt.ylabel('$\Delta$Ct')
plt.legend()

#plt.gcf()--current figure,plt.gca()currenr axis

#adjust the axis limits
#plt.xlim(-1,1)
plt.ylim(1,24)
#plt.axis('tight')


# In[4]:


#2D-scatter
data=pd.read_csv("mfa_2019_07_03_TRIPLEX_16_probe_concentration_A03_Amplitude.csv")
data=data.rename(columns={'Ch2 Amplitude':'ch2','Ch1 Amplitude':'ch1'})
data['label']=None
thresh1=[5000,13000,20500]
thresh2=[4000]
left=data[(data['ch2']<thresh2[0])]
right=data[(data['ch2']>=thresh2[0])]
neg=left[left['ch1']<thresh1[0]]
posb1=left[(left['ch1']<thresh1[1])&(left['ch1']>thresh1[0])]
posb2=left[(left['ch1']>thresh1[1])&(left['ch1']<thresh1[2])]
posb3=left[left['ch1']>thresh1[2]]
posg1=right[right['ch1']<thresh1[0]]
posg2=right[(right['ch1']>=thresh1[0])&(right['ch1']<thresh1[1])]
posg3=right[(right['ch1']>=thresh1[1])&(right['ch1']<thresh1[2])]
pos=right[right['ch1']>=thresh1[2]]
plt.figure(figsize=(12,8))
lst=[neg,posb1,posb2,posb3,posg1,posg2,posg3,pos]
label=['neg','pFAM1','pFAM2','pFAM3','pHEX','pH&F1','pH&F2','pH&F3']
clst=['gray','skyblue','dodgerblue','b','hotpink','violet','darkorchid','purple']
for i in range(len(lst)):
    plt.scatter(lst[i]['ch2'],lst[i]['ch1'],s=8,label=label[i],c=clst[i])
plt.legend(loc=2, prop={'size': 10},markerscale=4.5,fontsize='large')
plt.ylim(0,30000)


# In[5]:


#specificity of probes
#bar chart
fig=plt.figure(figsize=(12,5))#create bar figure
#plt.subplot(2,1,1)
obj=('Probe\nTarget','HPV16\nHPV16','HPV16\nHPV18','HPV16\nNeg','HPV18\nHPV16','HPV18\nHPV18','HPV18\nNeg')
dct=(0,20.37134,0,0,0,20.36165,0)
pos=np.arange(len(obj))
plt.bar(pos,dct,align='center',alpha=0.5)
plt.xticks(pos,obj)
plt.ylabel('$\Delta$Ct')
plt.title('Specificity of probes')
a=np.array(('HPV16','HPV16','HPV16','HPV18','HPV18','HPV18','HPV16','HPV18','Neg','HPV16','HPV18','Neg')).reshape((2,6))
df=pd.DataFrame(a,
             columns=['1','2','3','4','5','6'])
#plt.subplot(2,1,2)
#plt.table(cellText=df.values)


# In[6]:


#primer concentration optimization
data=pd.read_csv("ddPCR_evagreen_primer_concentration.csv")
primerc=data.iloc[:,0]
copyc=data.iloc[:,1]
resol=data.iloc[:,2]

plt.figure(figsize=(12,8))
f,ax1=plt.subplots()
color='steelblue'
ax1.set_xlabel('primer concentration(nM)')
ax1.set_ylabel('HPV16 concentration(copies/$\mu$L)',color=color)
ax1.set_ylim(0,1300)
#ax1.set_yticks(np.arange(-1,1001,20))
ax1.set_title('Primer Concentration Optimization (HPV16)')
ax1.tick_params(axis='y',labelcolor=color)

ax1.scatter(primerc,copyc,alpha=0.7,color=color)

color='firebrick'
ax2=ax1.twinx()
ax2.set_ylabel('Resolution',color=color)
#ax2.set_ylim(0,5)
ax2.plot(primerc,resol,color=color)
ax2.tick_params(axis='y',labelcolor=color)


# In[7]:


#primer concentration optimization
data=pd.read_csv("ddPCR_evagreen_primer_concentration.csv")
primerc=data.iloc[1:,0]
copyc=data.iloc[1:,3]
resol=data.iloc[1:,4]

plt.figure(figsize=(16,6))
f,ax1=plt.subplots()
color='steelblue'
ax1.set_xlabel('primer concentration(nM)')
ax1.set_ylabel('HPV18 concentration(copies/$\mu$L)',color=color)
ax1.set_ylim(0,1500)
#ax1.set_yticks(np.arange(-1,1001,20))
ax1.set_title('Primer Concentration Optimization(HPV18)')
ax1.tick_params(axis='y',labelcolor=color)

ax1.scatter(primerc,copyc,alpha=0.7,color=color)

color='firebrick'
ax2=ax1.twinx()
ax2.set_ylabel('Resolution',color=color)
ax2.set_ylim(0,6)
ax2.plot(primerc,resol,color=color)
ax2.tick_params(axis='y',labelcolor=color)


# In[8]:


#break axes
bax = brokenaxes(ylims=((-1, .6), (0.9, 1)), hspace=.05, despine=False)
x = np.linspace(0, 1, 100)
bax.plot(x, np.sin(10 * x), label='sin')
bax.plot(x, np.cos(10 * x), label='cos')
bax.legend(loc=3)
bax.set_xlabel('time')
bax.set_ylabel('value')

