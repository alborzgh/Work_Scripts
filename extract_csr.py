

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(1, './Utility_Scripts')

from deepsoil_output_plot import deepsoil_plotter
from kayen import liquefaction_triggering_Kayen

def extract_csr():
    directory = r"C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\PRDP\mean_profile_n0.15\Output\\"
    
    dir_list = os.listdir(directory)
    for cur_dir in dir_list:
        motions = os.listdir(directory + cur_dir + '/profile_0/')
        depth = np.zeros(0)
        results_list = []
        title = 'depth'
        for motion_ii, motion in enumerate(motions):
            print(directory + cur_dir + '/'+ motion + '/' + 'deepsoilout.db3')
            print("debug1")
            print(motion_ii)
            ds = deepsoil_plotter(directory + cur_dir + '/profile_0/'+ motion + '/' + 'deepsoilout.db3') 
            depth = ds.get_elem_mid_depth().squeeze()
            csr = ds.get_max_stress_ratio().squeeze()
            if motion_ii == 0:
                results_list.append(depth)
            title += ',' + motion
            results_list.append(csr.squeeze())
        to_write = np.stack(results_list, axis = 1)
        with open(directory+cur_dir + '.csv', 'w') as f:
            f.write(title + '\n')
        with open(directory+cur_dir + '.csv', 'ba') as f:
            np.savetxt(f, to_write, delimiter=',')

def decorate_plot(ax):
    ax.grid(which='both', lw=0.25, color=(0.25,0.25,0.25), linestyle=':')
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.invert_yaxis()

def print_profiles():
    directory = r"C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\PRDP\mean_profile_n0.15\Output\\"

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


        fig.savefig(directory+ cur_dir + '_profile.png')

def plot_profiles_fs():

    directory = r"C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\PRDP\mean_profile_n0.15\Output\\"

    profiles = [
        # 'Bhb-17-Simplified',
        # 'Bhd-01-Simplified',
        # 'Bhr-01-Simplified',
        # 'Bhr-04-Simplified',
        # 'Bhd-01-Simplified_new'
        'Bhr-04-Darendeli_soft2'
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
            csr = ds.get_max_stress_ratio().squeeze()*0.65
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



def plot_profiles_fs_excel():

    import openpyxl as xl
    directory = r"C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\PRDP\mean_profile_n0.15\Output\\"
    
    profiles = [
        # 'Bhb-17-Simplified_new',
        # # 'Bhd-01-Simplified',
        # 'Bhr-01-Simplified_new',
        # 'Bhr-04-Simplified_new',
        # 'Bhd-01-Simplified_new'
        'Darendeli',
        'Zhang'
    ]

    motions = [
        'Motion_Matched_Darfiel_N80E',
        'Motion_Matched_Duzce_531_N',
        'Motion_Matched_FKS_EW',
        'Motion_Matched_Hokkaid_63_BL',
        'Motion_Matched_Ibaraki_90',
        'Motion_matched_LANDERS_DSP090',
        'Motion_Matched_Manjil_L',
        'Motion_Matched_Romania_NE',
        'Motion_Matched_Sitka Alaska_90',
        'Motion_Matched_SOUTHERN_SUMATRA_N_BL',
    ]

    num_motions = len(motions)
    fig = plt.figure(figsize=(16,9), dpi=200)
    ax_pgd = fig.add_subplot(1,7,1)
    ax_pga = fig.add_subplot(1,7,2)
    ax_gmx = fig.add_subplot(1,7,3)
    ax_csr = fig.add_subplot(1,7,4)
    ax_fs  = fig.add_subplot(1,7,5)
    ax_pl  = fig.add_subplot(1,7,6)
    ax_vs  = fig.add_subplot(1,7,7)
    
    # wb = xl.load_workbook(r"C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\PRDP\mean_profile_n0.15\%s.xlsx"%profile, data_only=True)
    wb = xl.load_workbook(r"c:\Users\FEli\Desktop\PRDP\Bhr-04\Bhr-04-Darendeli_soft2.xlsx", data_only=True)

    sh = wb["Info"]
    waterlevel = float(sh[f'D7'].value)
    
    sh = wb["Deepsoil-Calibrated"]
    sublayer = 1
    Vs = []
    unit_weight = []
    in_depth = []
    pwp = []
    tot_stress = []
    eff_stress = []
    fc = []
    cur_depth = 0
    soil = []
    while True:

        if sh[f'A{sublayer + 1}'].value == None:
            break
        cur_depth = cur_depth + float(sh[f'C{sublayer + 1}'].value)/2.0
        
        in_depth.append(cur_depth)
        
        unit_weight.append(float(sh[f'E{sublayer + 1}'].value))			
        pwp.append(float(sh[f'R{sublayer + 1}'].value))
        tot_stress.append(float(sh[f'Q{sublayer + 1}'].value))
        eff_stress.append(float(sh[f'S{sublayer + 1}'].value))
        fc.append(float(sh[f'U{sublayer+1}'].value))
        Vs.append(float(sh[f'F{sublayer + 1}'].value))
        soil.append(str(sh[f'B{sublayer+1}'].value))
        
        cur_depth = cur_depth+float(sh[f'C{sublayer + 1}'].value)/2.0
        sublayer +=1
        
    tot_stress = np.array(tot_stress)
    eff_stress = np.array(eff_stress)
    Vs = np.array(Vs)
    fc = np.array(fc)
    csr_avg_all = 0
    pga_surface_avg = 0
    for profile in profiles:
        depth = np.zeros(0)
            
        for motion_ii, motion in enumerate(motions):

            # ds = deepsoil_plotter(os.path.join(directory,profile ,motion, 'deepsoilout.db3'))
            ds = deepsoil_plotter(os.path.join(directory,profile,'profile_0', motion, 'deepsoilout.db3'))
            depth = ds.get_elem_mid_depth().squeeze()
            pgd = ds.get_relPGD_profile().squeeze()
            pga = ds.get_totPGA_profile().squeeze()
            gmx = ds.get_max_strain().squeeze()
            csr = ds.get_max_stress_ratio().squeeze()*0.65
            
            pga_surface = pga[0]
            
            fs, pl, _, _ = liquefaction_triggering_Kayen(depth, Vs, csr=csr, sig_tot=tot_stress, sig_eff=eff_stress, Mw=8.2, fc=fc)
            
            # fs_simpl, pl_simpl, csr_simpl, crr_simpl = liquefaction_triggering_Kayen(depth, Vs, sig_tot=tot_stress, sig_eff=eff_stress, pga=0.24, Mw=8.2, fc=fc)
            fs_simpl, pl_simpl, csr_simpl, crr_simpl = liquefaction_triggering_Kayen(depth, Vs, sig_tot=tot_stress, sig_eff=eff_stress, pga=pga_surface, Mw=8.2, fc=fc)
            
            csr_avg_all += csr / num_motions/2
            pga_surface_avg += pga_surface/num_motions/2
            
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
			
            fs[depth<waterlevel]=np.NAN
            fs_simpl[depth<waterlevel]=np.NAN
            
            pl[depth<waterlevel]=np.NAN
            pl_simpl[depth<waterlevel]=np.NAN
            
            ax_fs.plot(fs, depth,'k-', lw=1.0, alpha=0.25)
            ax_fs.plot(fs_simpl, depth,'b-', lw=1.0, alpha=0.25)
            ax_pl.plot(pl, depth,'k-', lw=1.0, alpha=0.25)
            ax_pl.plot(pl_simpl, depth,'b-', lw=1.0, alpha=0.25)
            
        ax_pgd.plot(pgd_avg, depth, 'k-', lw=2.0)
        ax_pga.plot(pga_avg, depth, 'k-', lw=2.0)
        ax_gmx.plot(gmx_avg, depth, 'k-', lw=2.0)
        ax_csr.plot(csr_avg, depth, 'k-', lw=2.0)
        ax_csr.plot(csr_simpl_avg, depth, 'b-', lw=2.0)
        
        fs_avg[depth<waterlevel]=np.NAN
        fs_simpl_avg[depth<waterlevel]=np.NAN
        
        pl_avg[depth<waterlevel]=np.NAN
        pl_simpl_avg[depth<waterlevel]=np.NAN
        
        ax_fs.plot(fs_avg, depth,'k-', lw=2.0)
        ax_fs.plot(fs_simpl_avg, depth,'b-', lw=2.0)
        ax_pl.plot(pl_avg, depth,'k-', lw=2.0)
        ax_pl.plot(pl_simpl_avg, depth,'b-', lw=2.0)
        ax_vs.plot(Vs, depth, 'k-',lw=2.0)
        
        Vs_test = []
        ii = 0
        for idepth in depth:
            if(soil[ii].startswith('Compacted')):
                Vs_test.append(0.018*idepth**3.0 - 0.89*idepth**2.0 + 16.2*idepth + 104.0)
            elif(soil[ii].startswith('Stiff')):
                Vs_test.append(300)
            elif(soil[ii].startswith('Soft')):    
                Vs_test.append(1.5*idepth + 117.5)
            else:
                Vs_test.append(np.NAN)
            ii+=1
        ax_vs.plot(Vs_test, depth, 'g-',lw=2.0, label = "JD-2018")
        
        lines = [
                matplotlib.lines.Line2D([],[],color='k',linestyle='-', lw=2.0, label='Site-Specific'),
                matplotlib.lines.Line2D([],[],color='b',linestyle='-', lw=2.0, label='Simplified'),
                ]
        ax_pl.legend(handles=lines, loc='lower right')
        

        
    ## plot the average CSR  
    fs, pl, _, _ = liquefaction_triggering_Kayen(depth, Vs, csr=csr_avg_all, sig_tot=tot_stress, sig_eff=eff_stress, Mw=8.2, fc=fc)
    fs_simpl, pl_simpl, csr_simpl, crr_simpl = liquefaction_triggering_Kayen(depth, Vs, sig_tot=tot_stress, sig_eff=eff_stress, pga=pga_surface_avg, Mw=8.2, fc=fc)
 
    fs[depth<waterlevel]=np.NAN
    pl[depth<waterlevel]=np.NAN
        
    ax_fs.plot(fs, depth,'r-', lw=2.0)
    ax_pl.plot(pl, depth,'r-', lw=2.0)
    ax_fs.text(0, waterlevel-0.5,"Water Level")
    ax_fs.plot([0,2], [waterlevel,waterlevel],'k-.', lw=2.0)
    
    ax_pl.text(0, waterlevel-0.5,"Water Level")
    ax_pl.plot([0,1], [waterlevel,waterlevel],'k-.', lw=2.0)
    
    fs_simpl[depth<waterlevel]=np.NAN
    pl_simpl[depth<waterlevel]=np.NAN
    ax_fs.plot(fs_simpl, depth,'b--', lw=2.0)
    ax_pl.plot(pl_simpl, depth,'b--', lw=2.0)
    
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
    decorate_plot(ax_vs)
    ax_vs.set_xlabel("Vs")
    
    fig.suptitle(profile)
    fig.tight_layout(rect=[0.05,0.05,0.95,0.95])
        
    fig.savefig(os.path.join(directory, profile + '.png'))
    
    
		
if __name__ == "__main__":
    # extract_csr()
    # print_profiles()
    # plot_profiles_fs()
    plot_profiles_fs_excel()