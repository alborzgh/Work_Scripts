#TODO: - Fix the input generator
#      - Make functions respond to an array of element numbers or depths

import numpy as np
import matplotlib.pylab as plt
import sqlite3, matplotlib
import sys

from scipy.interpolate import interp1d

class deepsoil_plotter:
    def __init__(self, db_path=None):
        if db_path != None:
            self._ds_db_path = db_path
            self._read_db()

    def __del__(self):
        try:
            self._sql_connection.close()
        except:
            print("Something Happened!")

    def _read_db(self):
        if self._ds_db_path == None:
            print("No database is defined.")
            return
        try:
            self._sql_connection = sqlite3.connect(self._ds_db_path)
        except :
            print("There is an error here : " + str(sys.exc_info()[0]))
            raise
        self._sql_connection.row_factory = sqlite3.Row
        self._sql_cursor = self._sql_connection.cursor()

    def get_elem_number(self):
        self._sql_cursor.execute("SELECT LAYER_NUMBER from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_elem_top_depth(self):
        self._sql_cursor.execute("SELECT DEPTH_LAYER_TOP from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_elem_mid_depth(self):
        self._sql_cursor.execute("SELECT DEPTH_LAYER_MID from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_elem_bot_depth(self):
        bot_depth = 2.0*self.get_elem_mid_depth() - self.get_elem_top_depth()
        return bot_depth

    def get_node_depth(self):
        temp  = self.get_elem_top_depth()
        temp2 = self.get_elem_bot_depth()
        node_depth = np.append(temp, temp2[-1])
        return node_depth

    def get_totPGA_profile(self):
        self._sql_cursor.execute("SELECT PGA_TOTAL from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_relPGV_profile(self):
        self._sql_cursor.execute("SELECT PGV_RELATIVE from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_relminDisp_profile(self):
        self._sql_cursor.execute("SELECT MIN_DISP_RELATIVE from PROFILES")
        return np.array(self._sql_cursor.fetchall())
    
    def get_relmaxDisp_profile(self):
        self._sql_cursor.execute("SELECT MAX_DISP_RELATIVE from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_relPGD_profile(self):
        min_disp = np.abs(self.get_relminDisp_profile())
        max_disp = np.abs(self.get_relmaxDisp_profile())
        return np.maximum(min_disp, max_disp)

    def get_initial_eff_sig_v(self):
        self._sql_cursor.execute("SELECT INITIAL_EFFECTIVE_STRESS from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_max_strain(self):
        self._sql_cursor.execute("SELECT MAX_STRAIN from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_max_stress_ratio(self):
        self._sql_cursor.execute("SELECT MAX_STRESS_RATIO from PROFILES")
        return np.array(self._sql_cursor.fetchall())

    def get_input(self):
        self._sql_cursor.execute("SELECT INPUT from INPUT")
        return self._sql_cursor.fetchone()['INPUT'].rstrip(b'\0').decode("ascii")

    def generate_input_file(self, filename='deepsoil.dp'):
        # TODO: I need to add these lines before acceleration for the file to open:
        # [UNITS]:[ENGLISH]
        # [LAYER_NAMES]:[%num_Elems]
        #   [LAYER_NAMES]:[%elem_num][%elem_name]
        #   ...
        with open(filename, 'w', newline='\n') as f:
            f.write(self.get_input())
        return

    def get_FAS_frequencies(self):
        self._sql_cursor.execute("SELECT FREQUENCY from FOURIER_AMPLITUDE_SPECTRA")
        return np.array(self._sql_cursor.fetchall())

    def get_layer_FAS(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_FAS from FOURIER_AMPLITUDE_SPECTRA".format(layer_number))
        return np.array(self._sql_cursor.fetchall())

    def get_layer_FAS_ratio(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_FAS_RATIO from FOURIER_AMPLITUDE_SPECTRA".format(layer_number))
        return np.array(self._sql_cursor.fetchall())

    def get_input_FAS(self):
        self._sql_cursor.execute("SELECT INPUT_MOTION_FAS from FOURIER_AMPLITUDE_SPECTRA")
        return np.array(self._sql_cursor.fetchall())

    def get_response_spectrum_periods(self):
        self._sql_cursor.execute("SELECT PERIOD from RESPONSE_SPECTRA")
        return np.array(self._sql_cursor.fetchall())

    def get_response_spectrum(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_RS from RESPONSE_SPECTRA".format(layer_number))
        return np.array(self._sql_cursor.fetchall())

    def get_input_response_spectrum(self):
        self._sql_cursor.execute("SELECT INPUT_MOTION_RS from RESPONSE_SPECTRA")
        return np.array(self._sql_cursor.fetchall())

    def get_time_vector(self):
        self._sql_cursor.execute("SELECT TIME from TIME_HISTORIES")
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_acc_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_ACCEL from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_vel_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_VEL from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_disp_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_DISP from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_Arias_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_ARIAS from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_strain_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_STRAIN from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_stress_ratio_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_STRESS from TIME_HISTORIES".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_totVel_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_VEL_TOTAL from VEL_DISP".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_relVel_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_VEL_RELATIVE from VEL_DISP".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_totDisp_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_DISP_TOTAL from VEL_DISP".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_layer_relDisp_TH(self, layer_number=1):
        self._sql_cursor.execute("SELECT LAYER{0}_DISP_RELATIVE from VEL_DISP".format(layer_number))
        return np.array(self._sql_cursor.fetchall())
    
    def get_depth_response_spectrum(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_response_spectrum(layer_number=elemRecordNum) + portionOfBottom * self.get_response_spectrum(layer_number=elemRecordNum+1)
    
    def get_depth_acc_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_acc_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_acc_TH(layer_number=elemRecordNum+1)
    
    def get_depth_vel_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_vel_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_vel_TH(layer_number=elemRecordNum+1)
    
    def get_depth_disp_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_disp_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_disp_TH(layer_number=elemRecordNum+1)
    
    def get_depth_Arias_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_Arias_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_Arias_TH(layer_number=elemRecordNum+1)
    
    def get_depth_strain_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_strain_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_strain_TH(layer_number=elemRecordNum+1)
    
    def get_depth_stress_ratio_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_stress_ratio_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_stress_ratio_TH(layer_number=elemRecordNum+1)
    
    def get_depth_totVel_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_totVel_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_totVel_TH(layer_number=elemRecordNum+1)
    
    def get_depth_relVel_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_relVel_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_relVel_TH(layer_number=elemRecordNum+1)
    
    def get_depth_totDisp_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_totDisp_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_totDisp_TH(layer_number=elemRecordNum+1)
    
    def get_depth_relDisp_TH(self, depth=0.0):
        elem_top_depth = self.get_elem_top_depth().squeeze()
        elem_width = 2.0*(self.get_elem_mid_depth() - self.get_elem_top_depth())
        preElemNum = interp1d(elem_top_depth, self.get_elem_number().squeeze(), kind='previous', fill_value='extrapolate')
        elemRecordNum = int(preElemNum(depth))
        portionOfBottom = np.asscalar((depth - elem_top_depth[elemRecordNum-1])/elem_width[elemRecordNum-1])
        portionOfTop = 1 - portionOfBottom
        return portionOfTop * self.get_layer_relDisp_TH(layer_number=elemRecordNum) + portionOfBottom * self.get_layer_relDisp_TH(layer_number=elemRecordNum+1)
    
    
    
    
        

        

    



def test():
    test_database = r"./deepsoilout.db3"
    test_deepsoil_out = deepsoil_plotter(db_path=test_database)
    pgd_min = test_deepsoil_out.get_relminDisp_profile()
    pgd_max = test_deepsoil_out.get_relmaxDisp_profile()
    pgd = test_deepsoil_out.get_relPGD_profile()
    elem_num = test_deepsoil_out.get_elem_number()
    elem_top = test_deepsoil_out.get_elem_top_depth()
    elem_bot = test_deepsoil_out.get_elem_bot_depth()
    elem_mid = test_deepsoil_out.get_elem_mid_depth()
    
    time = test_deepsoil_out.get_time_vector()
    layer1_disp = test_deepsoil_out.get_layer_totDisp_TH(layer_number=1)
    layer10_disp = test_deepsoil_out.get_layer_totDisp_TH(layer_number=10)

    periods = test_deepsoil_out.get_response_spectrum_periods()
    input_rs = test_deepsoil_out.get_input_response_spectrum()
    layer1_rs = test_deepsoil_out.get_response_spectrum(layer_number=2)
    layer2_rs = test_deepsoil_out.get_response_spectrum(layer_number=3)

    depth_response_spectrum = test_deepsoil_out.get_depth_response_spectrum(depth=8.666666)

    # plt.figure()
    # plt.plot(pgd_min, elem_top)
    # plt.plot(pgd_max, elem_mid)
    # plt.plot(pgd, elem_bot)
    # plt.show()

    # plt.figure()
    # plt.plot(time, layer1_disp)
    # plt.plot(time, layer10_disp)
    # plt.show()

    plt.figure()
    # plt.semilogx(periods, input_rs)
    plt.semilogx(periods, layer1_rs)
    plt.semilogx(periods, layer2_rs)
    plt.semilogx(periods, depth_response_spectrum)
    plt.show()

    # with open('deepsoil.dp','w', newline='\n') as f:
    #     f.write(str(test_deepsoil_out.get_input()))

    # test_deepsoil_out.generate_input_file()


if __name__ == "__main__":
    test()
