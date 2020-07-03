# Comment by Feng

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from deepsoil_output_plot import deepsoil_plotter
from kayen import liquefaction_triggering_Kayen

def extract_csr():
    directory = r"C:\Golder-Projects\2020-PasayReclamation\Work Station\Output\\"
    dir_list = os.listdir(directory)
    for cur_dir in dir_list:
        motions = os.listdir(directory + cur_dir + '/profile_0/')
        depth = np.zeros(0)
        results_list = []
        title = 'depth'
        for motion_ii, motion in enumerate(motions):

            ds = deepsoil_plotter(directory + cur_dir + '/profile_0/' + motion + '/' + 'deepsoilout.db3')
            depth = ds.get_elem_mid_depth().squeeze()
            csr = ds.get_max_stress_ratio().squeeze()
            if motion_ii == 0:
                results_list.append(depth)
            title += ',' + motion
            results_list.append(csr.squeeze())
        to_write = np.stack(results_list, axis = 1)
        with open(cur_dir + '.csv', 'w') as f:
            f.write(title + '\n')
        with open(cur_dir + '.csv', 'ba') as f:
            np.savetxt(f, to_write, delimiter=',')

def decorate_plot(ax):
    ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.invert_yaxis()

def print_profiles():
    directory = r"C:\Golder-Projects\2020-PasayReclamation\Work Station\Output\\"
    dir_list = os.listdir(directory)
    for cur_dir in dir_list:
        motions = os.listdir(directory + cur_dir + '/profile_0/')
        depth = np.zeros(0)
        fig = plt.figure(figsize=(8,6), dpi=200, tight_layout=True)
        ax_pgd = fig.add_subplot(1,4,1)
        ax_pga = fig.add_subplot(1,4,2)
        ax_gmx = fig.add_subplot(1,4,3)
        ax_csr = fig.add_subplot(1,4,4)
        num_motions = len(motions)
        for motion_ii, motion in enumerate(motions):

            ds = deepsoil_plotter(directory + cur_dir + '/profile_0/' + motion + '/' + 'deepsoilout.db3')
            depth = ds.get_elem_mid_depth().squeeze()
            pgd = ds.get_relPGD_profile().squeeze()
            pga = ds.get_totPGA_profile().squeeze()
            gmx = ds.get_max_strain().squeeze()
            csr = ds.get_max_stress_ratio().squeeze()

            if motion_ii == 0:
                pgd_avg = pgd / num_motions
                pga_avg = pga / num_motions
                gmx_avg = gmx / num_motions
                csr_avg = csr / num_motions
            else:
                pgd_avg += pgd / num_motions
                pga_avg += pga / num_motions
                gmx_avg += gmx / num_motions
                csr_avg += csr / num_motions

            ax_pgd.plot(pgd, depth,'r-', lw=1.0, alpha=0.25)
            ax_pga.plot(pga, depth,'r-', lw=1.0, alpha=0.25)
            ax_gmx.plot(gmx, depth,'r-', lw=1.0, alpha=0.25)
            ax_csr.plot(csr, depth,'r-', lw=1.0, alpha=0.25)

        ax_pgd.plot(pgd_avg, depth, 'r-', lw=2.0)
        ax_pga.plot(pga_avg, depth, 'r-', lw=2.0)
        ax_gmx.plot(gmx_avg, depth, 'r-', lw=2.0)
        ax_csr.plot(csr_avg, depth, 'r-', lw=2.0)

        decorate_plot(ax_pgd)
        ax_pgd.set_xlabel("PGD (m)")
        ax_pgd.set_ylabel("Depth (m)")
        decorate_plot(ax_pga)
        ax_pga.set_xlabel("PGA (g)")
        decorate_plot(ax_gmx)
        ax_gmx.set_xlabel("$\gamma_{max}$ (%)")
        ax_gmx.set_xlim(right=1.0)
        decorate_plot(ax_csr)
        ax_csr.set_xlabel("CSR")


        fig.savefig(cur_dir + '.png')

def plot_profiles_fs():

    directory = r"C:\Golder-Projects\2020-PasayReclamation\Work Station\Output\\"
    
    profiles = [
        # 'Bhb-17-Simplified',
        # 'Bhd-01-Simplified',
        # 'Bhr-01-Simplified',
        # 'Bhr-04-Simplified',
        'Bhd-01-Simplified_new'
    ]

    motions = [
        'Matched_Darfiel_N80E',
        'Matched_Duzce_531_N',
        'Matched_FKS_EW',
        'Matched_Hokkaid_63_BL',
        'Matched_Ibaraki_90',
        'matched_LANDERS_DSP090',
        'Matched_Manjil_L',
        'Matched_Romania_NE',
        'Matched_Sitka Alaska_90',
        'Matched_SOUTHERN_SUMATRA_N_BL',
    ]

    num_motions = len(motions)

    for profile in profiles:
        depth = np.zeros(0)
        fig = plt.figure(figsize=(16,9), dpi=200)
        ax_pgd = fig.add_subplot(1,6,1)
        ax_pga = fig.add_subplot(1,6,2)
        ax_gmx = fig.add_subplot(1,6,3)
        ax_csr = fig.add_subplot(1,6,4)
        ax_fs  = fig.add_subplot(1,6,5)
        ax_pl  = fig.add_subplot(1,6,6)
        for motion_ii, motion in enumerate(motions):

            ds = deepsoil_plotter(os.path.join(directory,profile ,motion, 'deepsoilout.db3'))
            with open(os.path.join(directory,profile ,motion, 'deepsoilin.txt'), 'r') as f:
                ds_input = f.read().replace('\n',' ').replace('\t',' ').split()
            depth = ds.get_elem_mid_depth().squeeze()
            pgd = ds.get_relPGD_profile().squeeze()
            pga = ds.get_totPGA_profile().squeeze()
            gmx = ds.get_max_strain().squeeze()
            csr = ds.get_max_stress_ratio().squeeze()
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
                    
            fs, pl, _, _ = liquefaction_triggering_Kayen(depth, Vs, csr=0.65*csr, sig_tot=tot_stress, sig_eff=eff_stress, Mw=7.7, fc=0.0*Vs)
            fs_simpl, pl_simpl, csr_simpl, crr_simpl = liquefaction_triggering_Kayen(depth, Vs, sig_tot=tot_stress, sig_eff=eff_stress, pga=0.24, Mw=7.7, fc=0.0*Vs)

            if motion_ii == 0:
                pgd_avg = pgd / num_motions
                pga_avg = pga / num_motions
                gmx_avg = gmx / num_motions
                csr_avg = csr / num_motions
                fs_avg  = fs  / num_motions
                pl_avg  = pl  / num_motions
                csr_simpl_avg = csr_simpl / num_motions
                fs_simpl_avg  = fs_simpl  / num_motions
                pl_simpl_avg  = pl_simpl  / num_motions
            else:
                pgd_avg += pgd / num_motions
                pga_avg += pga / num_motions
                gmx_avg += gmx / num_motions
                csr_avg += csr / num_motions
                fs_avg  += fs  / num_motions
                pl_avg  += pl  / num_motions
                csr_simpl_avg += csr_simpl / num_motions
                fs_simpl_avg  += fs_simpl  / num_motions
                pl_simpl_avg  += pl_simpl  / num_motions

            ax_pgd.plot(pgd, depth,'k-', lw=1.0, alpha=0.25)
            ax_pga.plot(pga, depth,'k-', lw=1.0, alpha=0.25)
            ax_gmx.plot(gmx, depth,'k-', lw=1.0, alpha=0.25)
            ax_csr.plot(csr, depth,'k-', lw=1.0, alpha=0.25)
            ax_csr.plot(csr_simpl, depth,'b-', lw=1.0, alpha=0.25)
            ax_fs.plot(fs, depth,'k-', lw=1.0, alpha=0.25)
            ax_fs.plot(fs_simpl, depth,'b-', lw=1.0, alpha=0.25)
            ax_pl.plot(pl, depth,'k-', lw=1.0, alpha=0.25)
            ax_pl.plot(pl_simpl, depth,'b-', lw=1.0, alpha=0.25)

        ax_pgd.plot(pgd_avg, depth, 'k-', lw=2.0)
        ax_pga.plot(pga_avg, depth, 'k-', lw=2.0)
        ax_gmx.plot(gmx_avg, depth, 'k-', lw=2.0)
        ax_csr.plot(csr_avg, depth, 'k-', lw=2.0)
        ax_csr.plot(csr_simpl_avg, depth, 'b-', lw=2.0)
        ax_fs.plot(fs_avg, depth,'k-', lw=2.0)
        ax_fs.plot(fs_simpl_avg, depth,'b-', lw=2.0)
        ax_pl.plot(pl_avg, depth,'k-', lw=2.0)
        ax_pl.plot(pl_simpl_avg, depth,'b-', lw=2.0)

        decorate_plot(ax_pgd)
        ax_pgd.set_xlabel("PGD (m)")
        ax_pgd.set_ylabel("Depth (m)")
        decorate_plot(ax_pga)
        ax_pga.axvline(x=0.24, color='b', lw=2.0)
        ax_pga.set_xlabel("PGA (g)")
        decorate_plot(ax_gmx)
        ax_gmx.set_xlabel("$\gamma_{max}$ (%)")
        ax_gmx.set_xlim(right=1.0)
        decorate_plot(ax_csr)
        ax_csr.set_xlabel("CSR")
        decorate_plot(ax_fs)
        ax_fs.set_xlabel("FS")
        ax_fs.set_xlim(0.0,2.0)
        ax_fs.axvline(x=1.3, color='r', linestyle='--')
        decorate_plot(ax_pl)
        ax_pl.set_xlabel("PL")

        lines = [
                matplotlib.lines.Line2D([],[],color='k',linestyle='-', lw=2.0, label='Site-Specific'),
                matplotlib.lines.Line2D([],[],color='b',linestyle='-', lw=2.0, label='Simplified'),
                ]
        ax_pl.legend(handles=lines, loc='lower right')
        fig.suptitle(profile)
        fig.tight_layout(rect=[0.05,0.05,0.95,0.95])


        fig.savefig(os.path.join(directory, profile + '.png'))

if __name__ == "__main__":
    # extract_csr()
    # print_profiles()
    plot_profiles_fs()