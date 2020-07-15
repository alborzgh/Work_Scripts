import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(1, './Utility_Scripts')
from deepsoil_output_plot import deepsoil_plotter
from kayen import liquefaction_triggering_Kayen
from os_utils import ensure_dir

from multiprocessing import Pool
from scipy import interpolate
import glob  
import openpyxl as xl 


# TODO: Needs to be updatd for the dpz results

def run_parallel_aux(cwd):
    
    # print(os.path.join(cwd, 'deepsoilout.db3'))
    temp_file = os.path.join(cwd, 'deepsoilout.db3')
    if os.path.isfile(temp_file):
    
    
        ds = deepsoil_plotter(os.path.join(cwd, 'deepsoilout.db3'))
        depth = ds.get_elem_mid_depth().squeeze()
        pgd = ds.get_relPGD_profile().squeeze()
        pga = ds.get_totPGA_profile().squeeze()
        gmx = ds.get_max_strain().squeeze()
        csr = 0.65 * ds.get_max_stress_ratio().squeeze()
    
        with open(os.path.join(cwd, 'deepsoilin.txt'), 'r') as f:
            ds_input = f.read().replace('\n',' ').replace('\t',' ').split()
        Vs = []
        unit_weight = []
        in_depth = []
        cur_depth = 0.0
        layer_num = 0
        water_table = 1
        pwp = []
        tot_stress = []
        cur_pwp = 0.0
        cur_tot_stress = 0.0
        eff_stress = []
        for data in ds_input:
            if data.startswith('[WATER_TABLE]'):
                gwt = int(data.split(':')[1][1:-1])
            if data.startswith('[LAYER]'):
                if 'TOP' in data:
                    break
                layer_num = int(data.split(':')[1][1:-1])
            if data.startswith('[THICKNESS]'):
                thickness = float(data.split(':')[1][1:-1])
                in_depth.append(cur_depth + thickness / 2.0)
            if data.startswith('[WEIGHT]'):
                unit_weight.append(float(data.split(':')[1][1:-1]))
                cur_tot_stress += unit_weight[-1] * thickness / 2.0
                if layer_num > gwt:
                    cur_pwp += 9.81 * thickness / 2.0
                tot_stress.append(cur_tot_stress)
                pwp.append(cur_pwp)
                eff_stress.append(cur_tot_stress-cur_pwp)
                cur_tot_stress += unit_weight[-1] * thickness / 2.0
                if layer_num > gwt:
                    cur_pwp += 9.81 * thickness / 2.0
            if data.startswith('[SHEAR]'):
                Vs.append(float(data.split(':')[1][1:-1]))
        tot_stress = np.array(tot_stress)
        eff_stress = np.array(eff_stress)
        Vs = np.array(Vs)
                
        fs, pl, _, _ = liquefaction_triggering_Kayen(depth, Vs, csr=csr, sig_tot=tot_stress, sig_eff=eff_stress, Mw=7.7, fc=0.0*Vs)
        fs_simpl, pl_simpl, csr_simpl, crr_simpl = liquefaction_triggering_Kayen(depth, Vs, sig_tot=tot_stress, sig_eff=eff_stress, pga=0.24, Mw=7.7, fc=0.0*Vs)
    
        data_to_write = np.stack((depth, Vs, eff_stress, tot_stress, pgd, pga, gmx, csr, fs, pl, csr_simpl, crr_simpl, fs_simpl, pl_simpl), axis=1)
        with open(os.path.join(cwd, 'results.csv'), 'w') as f:
            f.write('depth, Vs, eff_stress, tot_stress, pgd, pga, gamma_max, csr, fs, pl, csr_simpl, crr_simpl, fs_simpl, pl_simpl\n')
        with open(os.path.join(cwd, 'results.csv'), 'ab') as f:
            np.savetxt(f, data_to_write, delimiter=',')
        del ds
        os.remove(os.path.join(cwd, 'deepsoilout.db3'))
        os.remove(os.path.join(cwd, 'deepsoilout_el.db3'))
        os.remove(os.path.join(cwd, 'strain.txt'))
        print(f'Done with {cwd}')

def run_parallel(folder_str):
    motions = [
        '../Motions/Matched_Darfiel_N80E.txt',
        '../Motions/Matched_Duzce_531_N.txt',
        '../Motions/Matched_FKS_EW.txt',
        '../Motions/Matched_Hokkaid_63_BL.txt',
        '../Motions/Matched_Ibaraki_90.txt',
        '../Motions/matched_LANDERS_DSP090.txt',
        '../Motions/Matched_Manjil_L.txt',
        '../Motions/Matched_Romania_NE.txt',
        '../Motions/Matched_Sitka Alaska_90.txt',
        '../Motions/Matched_SOUTHERN_SUMATRA_N_BL.txt',
    ]
    # folder_str = '../Bhr-04_parameter_study_temp2/Darendeli'
    cases = os.listdir(os.path.join(folder_str))
    commands = []
    for case in cases:
        if not case.startswith('profile'):
            continue
            
        for motion in motions:
            cwd = os.path.join(folder_str, case, motion.replace('../Motions/','Motion_').replace('.txt',''))
            commands.append(cwd)
    
    with Pool(8) as pool:
            pool.map(run_parallel_aux, commands)

def plot_parametric_study(profile_dir, list_ID_start, layer_depth):
    motions = [
        '../Motions/Matched_Darfiel_N80E.txt',
        '../Motions/Matched_Duzce_531_N.txt',
        '../Motions/Matched_FKS_EW.txt',
        '../Motions/Matched_Hokkaid_63_BL.txt',
        '../Motions/Matched_Ibaraki_90.txt',
        '../Motions/matched_LANDERS_DSP090.txt',
        '../Motions/Matched_Manjil_L.txt',
        '../Motions/Matched_Romania_NE.txt',
        '../Motions/Matched_Sitka Alaska_90.txt',
        '../Motions/Matched_SOUTHERN_SUMATRA_N_BL.txt',
    ]
    
	plots_dir = os.path.join(profile_dir, "FS") 
    ensure_dir(plots_dir)
	
    cases = os.listdir(os.path.join(profile_dir))
    # num_cases = len(cases)
    
    depth = np.zeros(0)
    Vs = np.zeros(0)
    cwd = os.path.join(profile_dir, 'Profile_1', motions[0].replace('../Motions/','Motion_').replace('.txt',''))
    # print(cwd)
    # print("debug0")
    res = np.loadtxt(os.path.join(cwd, 'results.csv'), delimiter=',', skiprows=1)
    indices = np.where(res[:,0] < layer_depth)
    cases_list = []
    num_cases = list_ID_start[1]-list_ID_start[0]
    
    ind = []
    for case in range(num_cases):
        # num_cases += 1
        
        this_case = {}
        # if not case.startswith('profile'):
        #     continue
        
        for ii in indices[0]:
            # print(ii)
            # print("ii")
            fs = np.zeros(0)
            pl = np.zeros(0)
            Vs = np.zeros(0)
            this_case[f'sublayer{ii}'] = {}
            ## add a new loop over the subfolders
            for ID_start in list_ID_start:                
                print({f'Profile_{ID_start+case}'})
                for motion in motions:
                    
                    # print(ID_start)
                    # print(case)
                    # print(num_cases)
                    cwd = os.path.join(profile_dir,  f'Profile_{ID_start+case}', motion.replace('../Motions/','Motion_').replace('.txt',''))
                    res = np.loadtxt(os.path.join(cwd, 'results.csv'), delimiter=',', skiprows=1)
                    # print(res[ii, 1])
                    Vs = np.append(Vs, res[ii, 1])
                    fs = np.append(fs, res[ii, 8])
                    pl = np.append(pl, res[ii, 9])
            # print(Vs)
            ## plot fs
            ser = pd.Series(fs)
            ser = ser.sort_values()
            # print(ser)
            count, division = pd.np.histogram(ser, density=1, bins=100)
            cdf = (count / count.sum()).cumsum()
            
            # print(division)
            if all(division <= 1.305):
                FS_13 = 0.0
            elif all(division > 1.295):
                FS_13 = 1.0
            else:                
                FS_13 = 1.0 - cdf[np.argmin(division <= 1.3)]
                
            if all(division <= 1.205):
                FS_12 = 0.0
            elif all(division > 1.195):
                FS_12 = 1.0
            else:                
                FS_12 = 1.0 - cdf[np.argmin(division <= 1.2)]
                
            
            if all(division <= 1.01):
                FS_10 = 0.0
            elif all(division > 0.99):
                FS_10 = 1.0
            else:                
                FS_10 = 1.0 - cdf[np.argmin(division <= 1.0)]
            # print(Vs[0])
            this_case[f'sublayer{ii}']['Vs'] = Vs[0]
            this_case[f'sublayer{ii}']['FS_13'] = FS_13
            this_case[f'sublayer{ii}']['FS_12'] = FS_12
            this_case[f'sublayer{ii}']['FS_10'] = FS_10
        
            fig = plt.figure(figsize=(5,4), dpi=200)
            ax = ser.hist(cumulative=True, density=1, bins=100)
            ax.set_xlim(0.0, 2.0)
            ax.set_ylim(0.0,1.0)
            ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
            ax.set_xlabel('FS')
            ax.set_ylabel('Cumulative Distribution')
            ax.set_title(f"case{case} for sublayer{ii+1}")
            
            ax.axvline(x=1.3, color='r', lw=1.0)
            ax.text(0.1, 0.9, f'Vs (m/s) = {round(Vs[0],0)}')
            ax.text(0.1, 0.8, f'Prob(FS>=1.3) = {round(FS_13,1)}')
            ax.text(0.1, 0.7, f'Prob(FS>=1.0) = {round(FS_10,1)}')
            fig.savefig(os.path.join(profile_dir,'FS', f'FS_distribution_layer{ii}_{case}.png'))
            # fig.savefig(os.path.join(profile_dir,'FS', f'FS_distribution_{case}_layer{ii}.png'))
            plt.close()
            
            ##plot lp
            PL_avg = sum(pl)/len(pl)
                
            this_case[f'sublayer{ii}']['PL_avg'] = PL_avg

        ind.append(ID_start)
        cases_list.append(this_case)


    for idx in indices[0]:
        # print(idx)
        Vs_list = []
        FS_13_list = []
        FS_12_list = []
        FS_10_list = []
        PL_avg_list = []
        for case in cases_list:
            # print(case[f'sublayer{idx}']['Vs'])
            Vs_list.append(case[f'sublayer{idx}']['Vs'])
            FS_13_list.append(case[f'sublayer{idx}']['FS_13'])
            FS_12_list.append(case[f'sublayer{idx}']['FS_12'])
            FS_10_list.append(case[f'sublayer{idx}']['FS_10'])
            PL_avg_list.append(case[f'sublayer{idx}']['PL_avg'])
        print(Vs_list)
        
        data_to_write = np.stack(( np.array(Vs_list), np.array(FS_13_list), np.array(FS_12_list), np.array(FS_10_list), np.array(PL_avg_list)), axis=1)
        
        with open(os.path.join(profile_dir, f'sublayer{idx}_VS_FS_PL.csv'), 'w') as f:
            f.write('Vs, FS_13, FS_12, FS_10, PL_avg\n')
        with open(os.path.join(profile_dir, f'sublayer{idx}_VS_FS_PL.csv'), 'ab') as f:
            np.savetxt(f, data_to_write, delimiter=',')
        
        
        
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)
    # ax.scatter(Vs, pl)
    # ax.set_ylim(0,1)
    # plt.show()

    # fig = plt.figure(figsize=(5,4), dpi=200)
    # ser = pd.Series(Vs)
    # ser = ser.sort_values()
    # ax = ser.hist(cumulative=True, density=1, bins=100)
    # ax.set_xlim(left = 0.0)
    # ax.set_ylim(0.0,1.0)
    # ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    # ax.set_xlabel('Vs (m/s)')
    # ax.set_ylabel('Cumulative Distribution')
    # fig.savefig(os.path.join(profile_dir, 'Vs_distribution.png'))
    

    # fig = plt.figure(figsize=(5,4), dpi=200)
    # ser = pd.Series(pl)
    # ser = ser.sort_values()
    # ax = ser.hist(cumulative=True, density=1, bins=100)
    # ax.set_xlim(0.0, 1.0)
    # ax.set_ylim(0.0,1.0)
    # ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    # ax.set_xlabel('PL')
    # ax.set_ylabel('Cumulative Distribution')
    # fig.savefig(os.path.join(profile_dir, 'PL_distribution.png'))


def plot_fs(profile_dir, waterlevel =[4.0]):
    
	
	plots_dir = os.path.join(profile_dir, "FS_VS") 
    ensure_dir(plots_dir)
	plots_dir = os.path.join(profile_dir, "PL") 
    ensure_dir(plots_dir)
	
    # sublayers = os.listdir(os.path.join(profile_dir))
    sublayers = glob.glob(f'{profile_dir}/*.csv')
    
    num_layers = len(sublayers)
    
    res = np.loadtxt(os.path.join(profile_dir,"profile_1\Motion_Matched_Darfiel_N80E", 'results.csv'), delimiter=',', skiprows=1)
    depth = res[:num_layers,0]
    
    # print(sublayers)
    FS_13_90th = []
    FS_13_50th = []
    
    FS_12_90th = []
    FS_10_50th = []    
    FS_10_84th = []
    FS_10_98th = []
    FS_10_90th = []
    
    PL_avg_all = []
    PL_avg15_all = []
    PL_avg0_all = []
    # num_layers = -1
    for isub in range(num_layers):
        print(f"%s  %s" %(num_layers,isub))
                
        # this_sub = {}
        # if not isub.startswith('sublayer'):
        #     continue
        
        # num_layers += 1
        res_unsorted = np.loadtxt(os.path.join(profile_dir, f'sublayer{isub}_VS_FS_PL.csv'), delimiter=',', skiprows=1)
        res = res_unsorted[np.argsort(res_unsorted[:, 0])]

        # res = np.loadtxt(isub, delimiter=',', skiprows=1)
                    # print(res[ii, 1])
        Vs_list = res[:, 0]
        FS_13_list = res[:, 1]
        FS_12_list = res[:, 2]
        FS_10_list = res[:, 3]
        PL_avg_list = res[:, 4]
        
        fig = plt.figure(figsize=(5,4), dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.plot(Vs_list, FS_13_list, 'k.',label="Limit = 1.3")
        # ax.plot(Vs_list, FS_12_list, 'b.',label="Limit = 1.2")
        ax.plot(Vs_list, FS_10_list, 'r-o',label="Limit = 1.0")
        legend = ax.legend(loc='lower right', shadow=False)#, fontsize='x-large')
        
        ax.set_xlabel('Vs (m/s)')
        ax.set_ylabel('P(FS>Limit)')
        # ax.set_xlim(left = 100.0)
        ax.set_title(f"sublayer{isub+1} at depth {round(depth[isub],0)}m")
        ax.set_xlim(100.0,300.0)
        ax.set_ylim(0.0,1.0)
        ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
        ax.axhline(y = 0.50, color='g', lw=0.25)
        # ax.axhline(y = 0.75, color='b', lw=0.25)
        ax.axhline(y = 0.90, color='g', lw=0.25)
        fig.savefig(os.path.join(profile_dir, "FS_VS",f'Vs_FS_Sublayer_{isub}.png'))
        plt.close()
        
        fig = plt.figure(figsize=(5,4), dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.plot(Vs_list, PL_avg_list, 'k.',label="average")
        # legend = ax.legend(loc='upper center', shadow=False, fontsize='x-large')
        
        ax.set_xlabel('Vs (m/s)')
        ax.set_ylabel('Average Probability of Liquefaction')
        # ax.set_xlim(left = 100.0)
        # ax.set_xlim(100.0,300.0)
        ax.set_ylim(0.0,1.0)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
        ax.axhline(y = 0.15, color='g', lw=0.25)
        # ax.axhline(y = 0.75, color='b', lw=0.25)
        # ax.axhline(y = 0.95, color='g', lw=0.25)
        ax.set_title(f"sublayer{isub+1} at depth {round(depth[isub],0)}m")
        fig.savefig(os.path.join(profile_dir,"PL", f'Vs_PL_Sublayer_{isub}.png'))
        plt.close()
        
        ### interpolate the FS plot
        poi = [0.5]
        # f = interpolate.interp1d(FS_13_list, Vs_list, fill_value=(-1), bounds_error=False)
        temp = interp([0.9], FS_13_list, Vs_list)
        FS_13_90th.append(temp[0])
        temp = interp([0.5], FS_13_list, Vs_list)
        FS_13_50th.append(temp[0])
        
        f = interpolate.interp1d(FS_12_list,Vs_list, fill_value=(-1), bounds_error=False)
        FS_12_90th.append(f([0.9])[0])
        
        
        f = interpolate.interp1d(FS_10_list, Vs_list, fill_value=(-1), bounds_error=False)
        FS_10_50th.append(f(poi)[0])
        FS_10_84th.append(f([0.84])[0])
        FS_10_98th.append(f([0.98])[0])  
        
        temp = interp([0.9], FS_10_list, Vs_list)
        FS_10_90th.append(temp[0])        
        
		f = interpolate.interp1d(PL_avg_list, Vs_list, fill_value=(-1), bounds_error=False)
        PL_avg_all.append(f(poi)[0])
        PL_avg15_all.append(f([0.15])[0])
        
        indices = np.where(PL_avg_list < 0.005)
        if len(indices[0])>0:
            PL_avg0_all.append(Vs_list[indices[0]][0])
        else:
            PL_avg0_all.append(-1)

    data_to_write = np.stack(( np.array(FS_13_90th), np.array(FS_12_90th), np.array(FS_10_50th),np.array(FS_10_84th),np.array(FS_10_98th), 
                              np.array(PL_avg_all),np.array(PL_avg15_all), np.array(PL_avg0_all)), axis=1)
    with open(os.path.join(profile_dir,"FS_VS", 'VS_FS_PL.csv'), 'w') as f:
        f.write(' FS_13, FS_12, FS_10_50th,FS_10_84th,FS_10_98th, PL_50, PL_15, PL_0\n')
    with open(os.path.join(profile_dir,"FS_VS", 'VS_FS_PL.csv'), 'ab') as f:
        np.savetxt(f, data_to_write, delimiter=',')    
        
    ## read idealized vs profile
    wb = xl.load_workbook(r'C:\Users\FEli\Desktop\PRDP\Bhr-04\Bhr-04-Darendeli_soft1.xlsx', data_only=True)

    # read idealized layering boundaries
    sh = wb['Deepsoil-Calibrated']
    layers_bottom = []
    Vs_list = []
    layer_bottom = 0.0
    row_num = 2
    while not sh[f'A{row_num}'].value is None:
       thickness = float(sh[f'C{row_num}'].value)
       layer_bottom += thickness
       layers_bottom.append(layer_bottom)
       Vs_list.append(sh[f'F{row_num}'].value)
       row_num += 1
        
    Vs_list = Vs_list[:num_layers]
    
    
    ### JD-2018 vs profile
    Vs_JD2018 = []
    for idepth in depth:
        Vs_test = 0.018*idepth**3.0 - 0.89*idepth**2.0 + 16.2*idepth + 104.0
        Vs_JD2018.append(Vs_test)
        
    fig = plt.figure(figsize=(5,4), dpi=200)
    ax = fig.add_subplot(1,1,1)
    FS_13_90th = np.array(FS_13_90th)
    FS_13_90th[depth<waterlevel]=np.NAN
    
    FS_13_50th = np.array(FS_13_50th)
    FS_13_50th[depth<waterlevel]=np.NAN
    
    FS_10_90th = np.array(FS_10_90th)
    FS_10_90th[depth<waterlevel]=np.NAN
    PL_avg15_all = np.array(PL_avg15_all)
    PL_avg15_all[depth<waterlevel]=np.NAN
    PL_avg0_all = np.array(PL_avg0_all)
    PL_avg0_all[depth<waterlevel]=np.NAN
    
    ax.plot(FS_13_90th, depth, 'b.',label="p(FS>=1.3)=90%")
    ax.plot(FS_13_50th, depth, 'b-',label="p(FS>=1.3)=50%")
    ax.plot(FS_10_90th, depth, 'g.',label="p(FS>=1.0)=90%")
    # ax.plot(PL_avg_all, depth, 'r-',label="PL=50%")
    ax.plot(PL_avg15_all, depth, 'r--',label="PL=15%")
    
    PL_avg0_all = np.array(PL_avg0_all)
    PL_avg0_all[PL_avg0_all < 0.0] = np.nan
    
    ax.plot(PL_avg0_all, depth, 'r-.',label="PL=0%")
    
    ax.plot(Vs_list, depth, 'k-',label="preliminary")
    # ax.plot(Vs_JD2018, depth, 'c-',label="JD2018")
    
    ax.set_ylim(ax.get_ylim()[::-1])
    
    legend = ax.legend(loc='bottom left', shadow=False)#, fontsize='x-large')
    
    ax.set_xlabel('Vs (m/s)')
    ax.set_ylabel('Depth (ft)')
    ax.set_title(f"Bhr-04")
    
    ax.set_xlim( 100.0,260)
    # ax.set_xlim(100.0,300.0)
    # ax.set_ylim(0.0,1.0)
    ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    # ax.axhline(y = 0.50, color='r', lw=0.25)
    # ax.axhline(y = 0.75, color='b', lw=0.25)
    # ax.axhline(y = 0.95, color='g', lw=0.25)
    ax.text(100,3.6,"water level", )
    ax.plot([100,300], [waterlevel,waterlevel], 'grey')
    fig.savefig(os.path.join(profile_dir, f'VS_profile.png'))
    plt.close()
    
    
    ### reliability plot
    fig = plt.figure(figsize=(5,4), dpi=200)
    ax = fig.add_subplot(1,1,1)
    
    ax.text(100,3.6,"water level", )
    ax.plot([100,300], [waterlevel,waterlevel], 'grey')
    
    FS_10_50th = np.array(FS_10_50th)
    FS_10_50th[depth<waterlevel]=np.NAN
    FS_10_84th = np.array(FS_10_84th)
    FS_10_84th[depth<waterlevel]=np.NAN
    FS_10_98th = np.array(FS_10_98th)
    FS_10_98th[depth<waterlevel]=np.NAN
    
    ax.plot(FS_10_50th, depth, 'g.',label="Median")
    ax.plot(FS_10_84th, depth, 'b.',label="Reliability index ~ 1.0") ## approximate the reliability index. needs confirm this.
    # ax.plot(FS_10_90th, depth, 'g.',label="Reliability index ~ 1.8??")
    ax.plot(FS_10_98th, depth, 'r.',label="Reliability index ~ 2.0")
    
    
    ax.plot(Vs_list, depth, 'k-',label="preliminary")
    # ax.plot(Vs_JD2018, depth, 'c-',label="JD2018")
    
    ax.set_ylim(ax.get_ylim()[::-1])
    
    legend = ax.legend(loc='bottom left', shadow=False)#, fontsize='x-large')
    
    ax.set_xlabel('Vs (m/s)')
    ax.set_ylabel('Depth (ft)')
    ax.set_title(f"Bhr-04")
    
    ax.set_xlim(100.0,260)
    # ax.set_xlim(100.0,300.0)
    # ax.set_ylim(0.0,1.0)
    ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    # ax.axhline(y = 0.50, color='r', lw=0.25)
    # ax.axhline(y = 0.75, color='b', lw=0.25)
    # ax.axhline(y = 0.95, color='g', lw=0.25)
    fig.savefig(os.path.join(profile_dir, f'VS_profile_reliability.png'))
    plt.close()

def interp(poi, FS_list, Vs_list):
    ##remove the ones with 0s and 1s FS
    indices = np.where(FS_list < 0.005)
    if(len(indices[0])):        
        ind_st = indices[0][len(indices[0])-1]
    else:
        ind_st = 0
    
    indices = np.where(FS_list > 0.995)
    if(len(indices[0])):        
        ind_ed = indices[0][0]+1
    else:
        ind_ed = len(FS_list)
        
    x = FS_list[ind_st:ind_ed]
    y = Vs_list[ind_st:ind_ed]
    f = interpolate.interp1d(x, y, fill_value=(-1), bounds_error=False)
    return(f(poi))
    
if __name__ == "__main__":

    # motions = [
    #     '../Motions/Matched_Darfiel_N80E.txt',
    #     '../Motions/Matched_Duzce_531_N.txt',
    #     '../Motions/Matched_FKS_EW.txt',
    #     '../Motions/Matched_Hokkaid_63_BL.txt',
    #     '../Motions/Matched_Ibaraki_90.txt',
    #     '../Motions/matched_LANDERS_DSP090.txt',
    #     '../Motions/Matched_Manjil_L.txt',
    #     '../Motions/Matched_Romania_NE.txt',
    #     '../Motions/Matched_Sitka Alaska_90.txt',
    #     '../Motions/Matched_SOUTHERN_SUMATRA_N_BL.txt',
    # ]
    # folder_str = '../Bhr-04_parameter_study_temp2/Darendeli'
    # cases = os.listdir(folder_str)
    
    # cwd = os.path.join(folder_str, cases[0], motions[1].replace('../Motions/','Motion_').replace('.txt',''))
    # print(cwd)
    # run_parallel_aux(cwd)
    
    folder_str = r'C:\Users\FEli\Documents\DEEPSOIL 7\Batch Output\Bhr04'
    
    list_ID_start = [1,41,81,121,161,201] 
    # run_parallel(folder_str)
    
    layer_depth = 15 #meter for Bhd-01
    layer_depth = 15 #meter for Bhr-04
    # plot_parametric_study(folder_str, list_ID_start, layer_depth=layer_depth)
    
    plot_fs(folder_str, waterlevel=[4.0])
    
    
    
    
    