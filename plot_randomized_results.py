import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from deepsoil_output_plot import deepsoil_plotter
from kayen import liquefaction_triggering_Kayen
from multiprocessing import Pool

# TODO: Needs to be updatd for the dpz results

def run_parallel_aux(cwd):
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

def run_parallel():
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

    cases = os.listdir('Bhr-01_parameter_study')
    commands = []
    for case in cases:
        if not case.startswith('Case'):
            continue
        for motion in motions:
            cwd = os.path.join('Bhr-01_parameter_study', case, motion.replace('../Motions/',''))
            commands.append(cwd)
    
    with Pool(8) as pool:
            pool.map(run_parallel_aux, commands)

def plot_parametric_study(profile_dir, layer_depth):
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

    cases = os.listdir(profile_dir)
    depth = np.zeros(0)
    Vs = np.zeros(0)
    cwd = os.path.join(profile_dir, 'Case1', motions[0].replace('../Motions/',''))
    res = np.loadtxt(os.path.join(cwd, 'results.csv'), delimiter=',', skiprows=1)
    indices = np.where(res[:,0] < layer_depth)
    cases_list = []
    num_cases = 0
    for case in cases:
        fs = np.zeros(0)
        Vs = np.zeros(0)
        this_case = {}
        if not case.startswith('Case'):
            continue
        for ii in indices[0]:
            num_cases += 1
            this_case[f'sublayer{ii}'] = {}
            for motion in motions:
                cwd = os.path.join(profile_dir, case, motion.replace('../Motions/',''))
                res = np.loadtxt(os.path.join(cwd, 'results.csv'), delimiter=',', skiprows=1)
                Vs = np.append(Vs, res[ii, 1])
                fs = np.append(fs, res[ii, 8])

            ser = pd.Series(fs)
            ser = ser.sort_values()
            count, division = pd.np.histogram(ser, density=1, bins=100)
            cdf = (count / count.sum()).cumsum()
            if all(division <= 1.3):
                FS_13 = 0.0
            elif all(division > 1.3):
                FS_13 = 1.0
            else:
                FS_13 = 1.0 - cdf[np.argmin(division <= 1.3)]

            this_case[f'sublayer{ii}']['Vs'] = Vs[0]
            this_case[f'sublayer{ii}']['FS_13'] = FS_13
        
            # fig = plt.figure(figsize=(5,4), dpi=200)
            # ax = ser.hist(cumulative=True, density=1, bins=100)
            # ax.set_xlim(0.0, 2.0)
            # ax.set_ylim(0.0,1.0)
            # ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
            # ax.set_xlabel('FS')
            # ax.set_ylabel('Cumulative Distribution')
            # ax.axvline(x=1.3, color='r', lw=1.0)
            # fig.savefig(os.path.join(profile_dir, case, f'FS_distribution_{case}_{ii}.png'))
            # plt.close()


        cases_list.append(this_case)

    for idx in indices[0]:
        Vs_list = []
        FS_13_list = []
        for case in cases_list:
            Vs_list.append(case[f'sublayer{idx}']['Vs'])
            FS_13_list.append(case[f'sublayer{idx}']['FS_13'])
        fig = plt.figure(figsize=(5,4), dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.plot(Vs_list, FS_13_list, 'k.')
        ax.set_xlabel('Vs (m/s)')
        ax.set_ylabel('P(FS>1.3)')
        ax.set_xlim(left = 0.0)
        ax.set_ylim(0.0,1.0)
        ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
        ax.axhline(y = 0.50, color='r', lw=0.25)
        ax.axhline(y = 0.75, color='b', lw=0.25)
        ax.axhline(y = 0.95, color='g', lw=0.25)
        fig.savefig(os.path.join(profile_dir, f'Vs_FS_Sublayer_{idx}.png'))
        plt.close()
            
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