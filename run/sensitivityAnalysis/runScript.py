import os

from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.BoundaryDict import BoundaryDict
import matplotlib.lines as mlines
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



myControl ={ 
'case1': {'E':2500000,'nu': 0.35 ,  'n' : 0.333 , 'k': 0.000001, 'relDensity' : 0.54 ,'k0' : 0.4, 'color': 'r','label': 'E=2.5*10^6' },
'case2': {'E':10000000,'nu': 0.35 ,  'n' : 0.333 , 'k': 0.000001, 'relDensity' : 0.54 ,'k0' : 0.4, 'color':'g', 'label': 'E=1*10^7'},
'case3': {'E':1000000,'nu': 0.35 ,  'n' : 0.333 , 'k': 0.000001, 'relDensity' : 0.54 ,'k0' : 0.4, 'color':'b','label': 'E=1*10^6' }

}
#case_list=['SumerExample']


fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()


d='SummerExample_b'
d1='SummerExample_p'

for i in myControl:
    file1=ParsedParameterFile(str(d+'/system/controlDict'))
    file1['endTime']= 13.7
    file1.writeFile()

    file2=ParsedParameterFile(str(d+'/constant/materialProperties'))
    file2['nu'][2]= myControl[i]['nu']
    file2['E'][2]= myControl[i]['E']
    file2['n'][2]= myControl[i]['n']
    file2['k'][2]= myControl[i]['k']  
    file2.writeFile()
    
    file4=ParsedParameterFile(str(d1+'/constant/materialProperties'))
    file4['relDensity'][2]= myControl[i]['relDensity']
    file4['k0'][2]= myControl[i]['k0']    
    file4.writeFile()
    print('BiotFoam starts')
    os.system(str('cd '+d+'; blockMesh > log.block; biotQSFoam > log.bio'));
    print('BiotFoam ends')
    os.system( 'cp -r ' +d+'/13.7/tauXZPrime2Mean '+d1+'/0/ ')
    os.system( 'cp -r ' +d+'/13.7/cv '+d1+'/0/ ')
    print('pressureBuildupFoam starts')
    file5=ParsedParameterFile(str(d1+'/system/controlDict'))
    file5['endTime']= 2500
    file5['writeFrequency']= 2500
    file5.writeFile()
    os.system(str('cd '+d1+'; blockMesh > log.block;  pressureBuildupFoam > log.pre'));
    print('pressureBuildupFoam ends')
    pE = pd.read_table(d1+'/postProcessing/Probes/0/pE', delim_whitespace=True, comment="#", names=['Time', 'z=0','z=-0.2','z=-0.4','z=-0.75','z=-1'])
    sigma0=np.array( [file4['gammaD'][2]*0.0001*(1+2*file4['k0'][2])/3, file4['gammaD'][2]*0.2*(1+2*file4['k0'][2])/3,file4['gammaD'][2]*0.4*(1+2*file4['k0'][2])/3,file4['gammaD'][2]*0.75*(1+2*file4['k0'][2])/3,file4['gammaD'][2]*1*(1+2*file4['k0'][2])/3])
   

    ax1.plot(pE['Time'],pE['z=0'],str(myControl[i]['color']+'-'))
    ax1.plot([0,2500],[sigma0[0],sigma0[0]],'k--' )
    ax2.plot(pE['Time'],pE['z=-0.2'],str(myControl[i]['color']+'-'))
    ax2.plot([0,2500],[sigma0[1],sigma0[1]],'k--' )
    ax3.plot(pE['Time'],pE['z=-0.4'],str(myControl[i]['color']+'-'))
    ax3.plot([0,2500],[sigma0[2],sigma0[2]],'k--' )
    ax4.plot(pE['Time'],pE['z=-0.75'],str(myControl[i]['color']+'-'))
    ax4.plot([0,2500],[sigma0[3],sigma0[3]],'k--' )
    ax5.plot(pE['Time'],pE['z=-1'],str(myControl[i]['color']+'-'))
    ax5.plot([0,2500],[sigma0[4],sigma0[4]],'k--' )
    

print([myControl[i]['label']  for i in myControl])
ax1.legend(handles=[mlines.Line2D([], [], color= myControl[i]['color'] , marker='.', linestyle='None',
    label=myControl[i]['label'])   for i in myControl],fontsize = 18)
ax2.legend(handles=[mlines.Line2D([], [], color= myControl[i]['color'] , marker='.', linestyle='None',
    label=myControl[i]['label'])   for i in myControl],fontsize = 18)
ax3.legend(handles=[mlines.Line2D([], [], color= myControl[i]['color'] , marker='.', linestyle='None',
    label=myControl[i]['label'])   for i in myControl],fontsize = 18)
ax4.legend(handles=[mlines.Line2D([], [], color= myControl[i]['color'] , marker='.', linestyle='None',
    label=myControl[i]['label'])   for i in myControl],fontsize = 18)
ax5.legend(handles=[mlines.Line2D([], [], color= myControl[i]['color'] , marker='.', linestyle='None',
    label=myControl[i]['label'])   for i in myControl],fontsize = 18)

ax1.set_ylabel('Accumulated pore pressure $[N/m^2]$')
ax1.set_xlabel('Time $[s]$') 
ax1.set_title('z= 0m',fontsize = 18) 
ax2.set_ylabel('Accumulated pore pressure $[N/m^2]$')
ax2.set_xlabel('Time $[s]$') 
ax2.set_title('z= 0.2m',fontsize = 18)
ax3.set_ylabel('Accumulated pore pressure $[N/m^2]$')
ax3.set_xlabel('Time $[s]$') 
ax3.set_title('z= 0.4m',fontsize = 18)
ax4.set_ylabel('Accumulated pore pressure $[N/m^2]$')
ax4.set_xlabel('Time $[s]$') 
ax4.set_title('z= 0.75m',fontsize = 18)
ax5.set_ylabel('Accumulated pore pressure $[N/m^2]$')
ax5.set_xlabel('Time $[s]$') 
ax5.set_title('z= 1m',fontsize = 18)

ax1.set_ylim(0, 20000)
ax2.set_ylim(0, 20000)
ax3.set_ylim(0, 20000)
ax4.set_ylim(0, 20000)
ax5.set_ylim(0, 20000)

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()

fig1.savefig('z0.png')
fig2.savefig('z02.png')
fig3.savefig('z04.png')
fig4.savefig('z075.png')
fig5.savefig('z1.png')
  
ax1.cla()
ax2.cla()
ax3.cla()
ax4.cla()
ax5.cla()

