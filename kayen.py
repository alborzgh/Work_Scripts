import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm

def liquefaction_triggering_Kayen(depth, Vs, Vs12=None, csr=None, fc=None, pga=None, Mw=None, sig_tot=None, sig_eff=None, K_sig=None, Patm=101.3, pl=0.15, dwf=None):
    Csv = (Patm/sig_eff)**0.25
    Csv = np.where(Csv>1.5,1.5, Csv)
    # print("CSV = ", Csv)
    Vs1 = Vs * Csv
    # print("Vs1 = ", Vs1)
    if csr is None:
        if pga is None:
            raise('Either CSR or PGA should be given!')
        if Vs12 is None:
            Vs12 = np.sum(np.where(depth <= 12.2, Vs/np.diff(np.append(depth, 12.2)), 0.0)) / np.sum(np.where(depth <= 12.2, np.diff(np.append(depth, 12.2)), 0.0))
        rd = (1 + (-23.013-2.949*pga+0.999*Mw+0.0525*Vs12)/(16.258 + 0.201*np.exp(0.341*(-depth+0.0785*Vs12+7.586))) ) / (1 + (-23.013-2.949*pga+0.999*Mw+0.0525*Vs12)/(16.258 + 0.201*np.exp(0.341*(0.0785*Vs12+7.586))))
        csr = 0.65 * pga * sig_tot / sig_eff * rd
    if dwf is None:
        dwf = 15*Mw**(-1.342)
    csr = csr/dwf
    # print("CSR = ", csr)
    pl_calc = norm.cdf(-((0.0073*Vs1)**2.8011 - 1.946*np.log(csr) - 2.6168*np.log(Mw) - 0.0099*np.log(sig_eff) + 0.0028*fc)/0.4809)
    crr = np.exp(((0.0073*Vs1)**2.8011 - 2.6168*np.log(Mw) - 0.0099*np.log(sig_eff) + 0.0028*fc - 0.4809 * norm.ppf(pl))/1.946)
    # print("CRR = ", crr)
    fs_liq = crr/csr

    return fs_liq, pl_calc, csr, crr

def kayen_test():
    import openpyxl as xl

    wb = xl.load_workbook(r"C:\Users\aghofrani\Downloads\Vs-Liq Calcs.xlsx", data_only=True)
    sh = wb["Kayen"]

    sublayer = 0
    depth = []
    csr = []
    Mw = []
    Vs = []
    fc = []
    Vs12 = []
    sig_tot = []
    sig_eff = []

    while True:
        if sh[f'A{sublayer + 6}'].value == None:
            break
        depth.append(float(sh[f'B{sublayer + 6}'].value))
        csr.append(float(sh[f'U{sublayer + 6}'].value))
        Mw.append(float(sh[f'D{sublayer + 6}'].value))
        Vs.append(float(sh[f'H{sublayer + 6}'].value))
        fc.append(float(sh[f'I{sublayer + 6}'].value))
        Vs12.append(float(sh[f'J{sublayer + 6}'].value))
        sig_tot.append(float(sh[f'L{sublayer + 6}'].value))
        sig_eff.append(float(sh[f'N{sublayer + 6}'].value))
        sublayer += 1
    
    depth = np.array(depth)
    csr = np.array(csr)
    Mw = np.array(Mw)
    Vs = np.array(Vs)
    fc = np.array(fc)
    Vs12 = np.array(Vs12)
    sig_tot = np.array(sig_tot)
    sig_eff = np.array(sig_eff)

    fs, _, _, _= liquefaction_triggering_Kayen(depth, Vs, Vs12, csr, fc, Mw=Mw, sig_tot=sig_tot, sig_eff=sig_eff, dwf=1.0)
    print("FS = ", fs)

if __name__ == "__main__":
    kayen_test()
    