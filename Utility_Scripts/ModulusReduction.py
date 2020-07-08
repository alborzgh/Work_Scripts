# Module for modulus reduction models
import numpy as np
from math import log, log10
from scipy import interpolate

def Darendeli(p, PI, OCR, frq=1, Ncyc=10, Patm=100.0):
    '''
    This function calculates the Darendeli modulus reduction and damping 
    curves for soil with a confining pressure of 'Conf_pressure', plastic 
    index of 'PI', overconsolidation ratio of 'OCR' and with a loading 
    frequency of 'frq' (default frq = 1) and number of cycles 'Ncyc' 
    (default Ncyc = 10), 'Patm' is used to normalize stresses (default Patm =
    100 kPa)
    Outputs are modulus reduction curve 'GGmax', damping ratio (%) 'Damp' and
    corresponding shear strain amplitudes 'gamma'
    
    written by : Pedro Arduino
    modified by: Alborz Ghofrani
    March 2015 - University of Washington

    revamped for python by Alborz Ghofrani
    November 2018
    '''

    numPoints = 31
    gamma   = np.logspace(-6,-1,numPoints)
    gamma_p = gamma * 100.0

    # constants
    f1  =  0.0352
    f2  =  0.0010
    f3  =  0.3246
    f4  =  0.3483
    f5  =  0.9190
    f6  =  0.8005
    f7  =  0.0129
    f8  = -0.1069
    f9  = -0.2889
    f10 =  0.2919
    f11 =  0.6329
    f12 = -0.0057

    a = f5
    b = f11 + f12 * log(Ncyc)
    cc1 = -1.1143*a*a+1.8618*a+0.2523
    cc2 =  0.0805*a*a-0.0710*a-0.0095
    cc3 = -0.0005*a*a+0.0002*a+0.0003

    # reference strain (%)
    gr = (f1 + f2 * PI * OCR ** f3) * (p / Patm) ** f4

    GGmax = np.zeros(numPoints)
    Damp  = np.zeros(numPoints)
    for ii in range(numPoints):
        # modulus ratio
        GGmax[ii] = 1.0 / (1.0 + (gamma_p[ii] / gr)**a)

        # damping ratio (%)
        Dmin = (f6 + f7 * PI * OCR**f8) * (p/Patm)**f9 * (1.0 + f10*log(frq))
        DMa1 = 100.0 * ( 4.0*( gamma_p[ii] - gr * log( (gamma_p[ii]+gr) / gr ) ) / \
                        (gamma_p[ii]**2.0 / ( gamma_p[ii]+gr) ) - 2.0 ) / np.pi
        Dmasing = cc1 * DMa1 + cc2 * DMa1**2.0 + cc3 * DMa1**3.0

        Damp[ii] = b * (GGmax[ii]**0.1) * Dmasing + Dmin

    return GGmax, Damp, gamma


def Darendeli_Corrected(p, PI, OCR, Gmax, tau_f, frq=1, Ncyc=10, Patm=100.0, gamma_1=-1, gamma_2=0.08):
    '''
    This function calculates the Darendeli modulus reduction and damping 
    curves for soil with a confining pressure of 'Conf_pressure', plastic 
    index of 'PI', overconsolidation ratio of 'OCR' and with a loading 
    frequency of 'frq' (default frq = 1) and number of cycles 'Ncyc' 
    (default Ncyc = 10), 'Patm' is used to normalize stresses (default Patm =
    100 kPa)
    This function also applies correction to the Darendeli curves to account
    for inconsistencies with the strength (tau_f) of the material. Initial
    stiffness of the material (Gmax) and a transition strain (gamma_1) and a
    strain at which 0.95 tau_f is enforced (gamma_2, default = 0.08), is 
    required to apply the correction. The correction uses a second hyperbola
    starting at gamma_1 to enforce 0.95 tau_f at gamma_2.
    Outputs are modulus reduction curve 'GGmax', damping ratio (%) 'Damp' and
    corresponding shear strain amplitudes 'gamma'
    
    written by : Pedro Arduino
    modified by: Alborz Ghofrani
    March 2015 - University of Washington

    revamped for python by Alborz Ghofrani
    November 2018
    '''

    numPoints = 31
    gamma   = np.logspace(-6,-1,numPoints)
    gamma_p = gamma * 100.0

    # constants
    f1  =  0.0352
    f2  =  0.0010
    f3  =  0.3246
    f4  =  0.3483
    f5  =  0.9190
    f6  =  0.8005
    f7  =  0.0129
    f8  = -0.1069
    f9  = -0.2889
    f10 =  0.2919
    f11 =  0.6329
    f12 = -0.0057

    a = f5
    b = f11 + f12 * log(Ncyc)
    cc1 = -1.1143*a*a+1.8618*a+0.2523
    cc2 =  0.0805*a*a-0.0710*a-0.0095
    cc3 = -0.0005*a*a+0.0002*a+0.0003

    # reference strain (%)
    gr = (f1 + f2 * PI * OCR ** f3) * (p / Patm) ** f4

    GGmax = np.zeros(numPoints)
    Damp  = np.zeros(numPoints)
    for ii in range(numPoints):
        # modulus ratio
        GGmax[ii] = 1.0 / (1.0 + (gamma_p[ii] / gr)**a)

        # damping ratio (%)
        Dmin = (f6 + f7 * PI * OCR**f8) * (p/Patm)**f9 * (1.0 + f10*log(frq))
        DMa1 = 100.0 * ( 4.0*( gamma_p[ii] - gr * log( (gamma_p[ii]+gr) / gr ) ) / \
                        (gamma_p[ii]**2.0 / ( gamma_p[ii]+gr) ) - 2.0 ) / np.pi
        Dmasing = cc1 * DMa1 + cc2 * DMa1**2.0 + cc3 * DMa1**3.0

        Damp[ii] = b * (GGmax[ii]**0.1) * Dmasing + Dmin

    GGmax_org = np.array(GGmax)
    
    if (gamma_2 < 0):
        gamma_2 = 0.08

    if (gamma_1 < 0):
        tau_max = GGmax[-1] * Gmax * gamma[-1]
        tau_1   = 0.5 * tau_max
        if tau_1 > tau_f:
            tau_1 = 0.7 * tau_f
        
        gamma_1 = interpolate.pchip_interpolate(GGmax*Gmax*gamma, gamma, tau_1)

    else:
        tau_1 = interpolate.pchip_interpolate(gamma, GGmax*Gmax*gamma, gamma_1)

    tau_m = 0.99 * tau_f - tau_1
    gamma_m = gamma_2

    Gmax_p = Gmax * (1+(1-a)*(gamma_1*100/gr)**a) / (1.0+(gamma_1*100/gr)**a)**2
    gr_p   = gamma_1 + tau_m*(gamma_1-gamma_m) / (tau_m+Gmax_p*(gamma_1-gamma_m))

    GGmax[np.where(gamma>gamma_1)] = ( tau_1 + Gmax_p*( gamma[ np.where(gamma>gamma_1) ] - gamma_1 ) / \
        ( 1.0 + ( gamma[ np.where(gamma>gamma_1) ]-gamma_1 ) / ( gr_p-gamma_1 ) ) ) / \
        Gmax / gamma[np.where(gamma>gamma_1)]


    return GGmax, Damp, gamma, GGmax_org


def Menq(p, Cu, D50, Ncyc=10, Patm=100.0):
    '''
    This function calculates the Menq modulus reduction and damping 
    curves for soil with a confining pressure of 'Conf_pressure', uniformuty 
    coefficient of 'Cu', 'D50 (mm)' and number of cycles 'Ncyc' (default 
    Ncyc = 10), 'Patm' is used to normalize stresses (default Patm = 100 kPa)
    Outputs are modulus reduction curve 'GGmax', damping ratio (%) 'Damp' and
    corresponding shear strain amplitudes 'gamma'
    
    written by: Alborz Ghofrani
    Dec 2015 - University of Washington

    revamped for python by Alborz Ghofrani
    November 2018
    '''
    
    numPoints = 31
    gamma   = np.logspace(-6,-1,numPoints)
    gamma_p = gamma * 100.0
    
    # constants
    f1  =  0.1200
    f2  = -0.6000
    f3  =  0.5000
    f4  = -0.1500
    f5  =  0.8600
    f6  =  0.1000
    f7  =  0.5500
    f8  =  0.1000
    f9  = -0.3000
    f10 = -0.0800
    f11 =  0.6329
    f12 = -0.0057

    # 10 based for a, natural based for damping parameters!!
    a   = f5+f6*log10(p/Patm)
    b   = f11+f12*log(Ncyc)
    cc1 = -1.1143*a*a+1.8618*a+0.2523
    cc2 =  0.0805*a*a-0.0710*a-0.0095
    cc3 = -0.0005*a*a+0.0002*a+0.0003

    # reference strain (%)
    gr = f1*Cu**(f2)*(p/Patm)**(f3*Cu**(f4))

    GGmax = np.zeros(numPoints)
    Damp  = np.zeros(numPoints)
    for ii in range(numPoints):
        # modulus ratio 
        GGmax[ii] = 1.0/(1.0+(gamma_p[ii]/gr)**a)

        # damping ratio (%)
        Dmin = f7*Cu**f8 * D50**f9 * (p/Patm)**(f10)
        DMa1 = 100.0*(4.0*(gamma_p[ii] - gr*log((gamma_p[ii]+gr)/gr))/ \
                        (gamma_p[ii]**2.0 / (gamma_p[ii]+gr))-2.0)/np.pi
        Dmasing = cc1*DMa1+cc2*DMa1**2.0+cc3*DMa1**3.0
        Damp[ii] = b * (GGmax[ii]**0.1) * Dmasing + Dmin 
    
    return GGmax, Damp, gamma


def Menq_Corrected(p, Cu, D50, Gmax, tau_f, Ncyc=10, Patm=100, gamma_1=-1, gamma_2=0.08):
    '''
    This function calculates the Menq modulus reduction and damping 
    curves for soil with a confining pressure of 'Conf_pressure', uniformuty 
    coefficient of 'Cu', 'D50 (mm)' and number of cycles 'Ncyc' (default 
    Ncyc = 10), 'Patm' is used to normalize stresses (default Patm = 100 kPa)
    This function also applies correction to the Darendeli curves to account
    for inconsistencies with the strength (tau_f) of the material. Initial
    stiffness of the material (Gmax) and a transition strain (gamma_1) and a
    strain at which 0.95 tau_f is enforced (gamma_2, default = 0.08), is 
    required to apply the correction. The correction uses a second hyperbola
    starting at gamma_1 to enforce 0.95 tau_f at gamma_2.
    Outputs are modulus reduction curve 'GGmax', damping ratio (%) 'Damp' and
    corresponding shear strain amplitudes 'gamma'   
    
    written by: Alborz Ghofrani
    Dec 2015 - University of Washington

    revamped for python by Alborz Ghofrani
    November 2018
    '''
    
    numPoints = 31
    gamma   = np.logspace(-6,-1,numPoints)
    gamma_p = gamma * 100.0
    
    # constants
    f1  =  0.1200
    f2  = -0.6000
    f3  =  0.5000
    f4  = -0.1500
    f5  =  0.8600
    f6  =  0.1000
    f7  =  0.5500
    f8  =  0.1000
    f9  = -0.3000
    f10 = -0.0800
    f11 =  0.6329
    f12 = -0.0057

    # 10 based for a, natural based for damping parameters!!
    a   = f5+f6*log10(p/Patm)
    b   = f11+f12*log(Ncyc)
    cc1 = -1.1143*a*a+1.8618*a+0.2523
    cc2 =  0.0805*a*a-0.0710*a-0.0095
    cc3 = -0.0005*a*a+0.0002*a+0.0003

    # reference strain (%)
    gr = f1*Cu**(f2)*(p/Patm)**(f3*Cu**(f4))

    GGmax = np.zeros(numPoints)
    Damp  = np.zeros(numPoints)
    for ii in range(numPoints):
        # modulus ratio 
        GGmax[ii] = 1.0/(1.0+(gamma_p[ii]/gr)**a)

        # damping ratio (%)
        Dmin = f7*Cu**f8 * D50**f9 * (p/Patm)**(f10)
        DMa1 = 100.0*(4.0*(gamma_p[ii] - gr*log((gamma_p[ii]+gr)/gr))/ \
                        (gamma_p[ii]**2.0 / (gamma_p[ii]+gr))-2.0)/np.pi
        Dmasing = cc1*DMa1+cc2*DMa1**2.0+cc3*DMa1**3.0
        Damp[ii] = b * (GGmax[ii]**0.1) * Dmasing + Dmin 
    
    GGmax_org = np.array(GGmax)

    if gamma_2 < 0:
        gamma_2 = 0.08
    
    if gamma_1 < 0:
        tau_max = GGmax[-1] * Gmax * gamma[-1]
        tau_1   = 0.5 * tau_max
        if (tau_1 > tau_f):
            tau_1 = 0.7 * tau_f

        gamma_1 = interpolate.pchip_interpolate(GGmax*Gmax*gamma,gamma, tau_1)

    else:
        tau_1   = interpolate.pchip_interpolate(gamma,GGmax*Gmax*gamma, gamma_1)

    tau_m = 0.99 * tau_f - tau_1
    gamma_m = gamma_2

    Gmax_p = Gmax * (1.0+(1.0-a)*(gamma_1*100.0/gr)**a) / (1.0+(gamma_1*100.0/gr)**a)**2.0
    gr_p   = gamma_1 + tau_m*(gamma_1-gamma_m) / (tau_m+Gmax_p*(gamma_1-gamma_m))

    GGmax[ np.where(gamma>gamma_1) ] = (tau_1 + Gmax_p*(gamma[ np.where(gamma>gamma_1) ] - gamma_1) / \
        (1.0+(gamma[ np.where(gamma>gamma_1) ] - gamma_1) / (gr_p-gamma_1))) / \
        Gmax / gamma[ np.where(gamma>gamma_1)]

    return GGmax, Damp, gamma, GGmax_org

def VuceticDobry(PI):
    '''
    This function returns the Vucetic-Dobry modulus reduction and damping 
    curves for soil with a plasticity index of PI   
    
    written by: Alborz Ghofrani
    November 2018
    '''

    vd_PI = np.array([0,15,30,50,100,200])
    
    vd_MR_gamma = np.array([ \
            [0.0001,0.000134861191225868,0.000181875408988603,0.000245279343108951,0.000330786643947717,0.000446102808443964,0.000601619561559586,0.000811351307367107,0.00109419803814194,0.0014756485086083,0.00199007715701592,0.00268384176026558,0.00361946096851148,0.00488124817808995,0.00658290943966313,0.00887779008764986,0.0119726934667367,0.0161465170310628,0.0217753852095791,0.0293665438876601,0.0396040709087659,0.0534105218014995,0.0720300659414546,0.0971406049694244,0.131004977025781,0.176674872582144,0.238265837761065,0.321328147088868,0.433346966908059,0.58441688171338,0.788151568403742],\
            [0.0001,0.000136130192667436,0.000185314293556734,0.000252268704859083,0.000343413873964317,0.000467489968274334,0.000636394994512791,0.000866325732156187,0.0011793308883118,0.00160542541044545,0.00218546870437138,0.0029750827579473,0.00404998589040935,0.00551325359561826,0.00750520274195955,0.0102168469527112,0.013908213441263,0.0189332777541897,0.0257739074850394,0.0350860699173111,0.0477627345778673,0.0650195026040873,0.0885111741663532,0.120490431924867,0.164023857125148,0.223285992725026,0.303959652095977,0.413780860029525,0.563280681979171,0.766795077636698,1.04383961655126],\
            [0.0001,0.000136859308359272,0.000187304702845783,0.000256343920839128,0.000350830517081471,0.00048014421919096,0.000657122057511775,0.000899332702986834,0.00123082051715652,0.00168449244692444,0.00230538471222496,0.00315513357217148,0.00431809398468509,0.00590971356174334,0.00808799310661606,0.0110691714258603,0.0151491914545346,0.020733078646698,0.0283751480374548,0.0388340313499803,0.0531479867136059,0.072737967023119,0.0995486785824361,0.136241632988716,0.186459356605735,0.255186985821758,0.349247143818532,0.477977225494556,0.654156324926688,0.895273821883099,1.22526556055083],\
            [0.0001,0.000137592329207499,0.000189316490567449,0.00026048496894565,0.000358407336007752,0.000493141001663615,0.000678524190466162,0.000933597237898723,0.00128455818504173,0.0017674535268245,0.00243188047521793,0.00334608098939477,0.00460395077047761,0.00633468310066677,0.00871603802612126,0.0119925997347516,0.0165008973075771,0.022703968945633,0.0312389196948438,0.042982357227396,0.0591404264574622,0.0813726902660706,0.11196257987589,0.154051921492045,0.211963626969814,0.291645691420463,0.401282099858732,0.552133387888393,0.759693188727919,1.04527955320147,1.43822448397965],\
            [0.0001,0.000138493245945231,0.000191803791724465,0.000265635297005243,0.000367886945198819,0.000509498571814602,0.000705621110150641,0.000977237579522402,0.00135340804447719,0.00187437873168035,0.00259588794681119,0.00359512947863985,0.00497901151090223,0.00689559465743523,0.00954993286830805,0.0132260120149104,0.0183171333485558,0.0253679925385315,0.0351329562977565,0.0486567715732828,0.0673863423239963,0.0933255328082681,0.129249559681852,0.179001910573317,0.247905556256969,0.343332451738859,0.475492256796493,0.658524660655702,0.912012177891905,1.2630752685783,1.74927393818555],\
            [0.0001,0.000139151071149167,0.000193630206019606,0.000269438505744622,0.000374926566831953,0.000521714333769463,0.000725971083778952,0.00101019653931163,0.00140569930516395,0.00195604564027205,0.00272185846060515,0.00378749520309631,0.00527034014483186,0.00733373476473811,0.0102049704803719,0.0142003257338939,0.0197599053653843,0.0274961199739941,0.038261145468273,0.053240793753043,0.0740851347956785,0.103090258630491,0.143451199134775,0.199613880172366,0.277764852422263,0.386512767421483,0.537836655995285,0.748405467850302,1.04141422505263,1.44913904926054,2.01649250948691],\
        ])
    vd_MR_GGmax = np.array([\
            [1,0.99999999406266,0.99999999359211,0.999055528388734,0.996712317389739,0.993023023750591,0.985412496447743,0.975148848668553,0.962009685740508,0.944649916547431,0.921238325019586,0.897949418708406,0.868345019491325,0.828420174051907,0.789209711367892,0.73883176717315,0.692472187375542,0.63240778292714,0.571667761859095,0.503870865958534,0.440820779803125,0.379049601594102,0.317168357440657,0.26925261902455,0.222477202027848,0.182223044395537,0.140285201365475,0.108885957810584,0.081458583138562,0.0622865624321664,0.0372703656447886],\
            [1,1,1,1,1,1,1,0.997889711558759,0.989493803592103,0.978387114833688,0.963483938932076,0.944459020615652,0.920033736475963,0.889277793910697,0.859738903656562,0.819527556532453,0.784770792394652,0.737456342408898,0.688086252102765,0.635104889857265,0.57907444262133,0.512414085296533,0.449014819320199,0.384472786808132,0.328768992865302,0.278446542279625,0.230128052064018,0.189928003792894,0.151512958252934,0.122784132156016,0.0945931984538658],\
            [1,1,1,1,1,1,1,1,1,0.997517462860287,0.992008229422519,0.982780911982445,0.969946588004196,0.952381620974976,0.928342328585252,0.899847227229121,0.862966591988467,0.825290696520132,0.785758430476676,0.731654844590854,0.664517852511406,0.617718510949757,0.553669255771791,0.499714159617923,0.437637835895207,0.385513267799707,0.332953896163892,0.281363917593255,0.233145673369166,0.189232236421708,0.143307107654254],\
            [1,1,1,1,1,1,1,1,1,1,0.998896105865268,0.995007975223537,0.989093286246614,0.979332745077274,0.966096054508625,0.948989367353943,0.925830106783811,0.897600740327169,0.861382470279891,0.822889723460774,0.769926656536405,0.720044880664534,0.668299268782147,0.616217949753819,0.552965121047496,0.496123073893052,0.433857987753117,0.375566942872395,0.311321553720572,0.250237220970657,0.186246738824423],\
            [1,0.999999999934021,0.999999999842645,0.999999999716096,0.999999999540834,0.999999999298107,0.999999998961948,0.999999998496389,0.999999997851622,0.999999996958664,0.999999995721976,0.999999994009247,0.99999999359211,0.99999999359211,0.991333039268486,0.979280623466247,0.964766905183548,0.950162088628382,0.929935404116786,0.901922812189575,0.872087847565714,0.833971725526782,0.792154590141348,0.741701154263442,0.681702556683888,0.61655553134187,0.550361426162161,0.477104499019744,0.399266490244841,0.322484804449872,0.249238862790442],\
            [1,0.999999999901313,0.999999999763989,0.999999999572902,0.999999999307002,0.999999998937,0.999999998422137,0.9999999977057,0.999999996708771,0.999999995321533,0.99999999359211,0.99999999359211,0.99999999359211,0.999022190223907,0.996528298492521,0.993058021434996,0.988058837059477,0.977844271794707,0.964077090001611,0.948865487323859,0.927698379258811,0.901687537975985,0.86649031891685,0.824590173831254,0.772477840265678,0.712010598034234,0.640189543751187,0.565238186607509,0.47275481378168,0.379371788262037,0.292283480833888],\
        ])

    vd_damp_gamma = np.array([\
            [0.0001,0.000137936822210558,0.000190265669215471,0.000262446417873473,0.000362010248820111,0.000499345433298996,0.000688781222546178,0.00095008293036323,0.00131051420250798,0.00180768164555755,0.00249345861756562,0.00343939758020533,0.00474419572532207,0.0065439928229584,0.00902657574567582,0.0124509717380142,0.0171744747497515,0.0236899247011619,0.0326771293168568,0.045073793769307,0.0621733587751226,0.0857599553559733,0.118294557147223,0.163171752976932,0.225073930801642,0.310459827772175,0.42823842066931,0.590698468955929,0.814790696924229,1.12389639500454,1.55026697220829],\
            [0.0001,0.000137936822210558,0.000190265669215471,0.000262446417873473,0.000362010248820111,0.000499345433298996,0.000688781222546178,0.00095008293036323,0.00131051420250798,0.00180768164555755,0.00249345861756562,0.00343939758020533,0.00474419572532207,0.0065439928229584,0.00902657574567582,0.0124509717380142,0.0171744747497515,0.0236899247011619,0.0326771293168568,0.045073793769307,0.0621733587751226,0.0857599553559733,0.118294557147223,0.163171752976932,0.225073930801642,0.310459827772175,0.42823842066931,0.590698468955929,0.814790696924229,1.12389639500454,1.55026697220829],\
            [0.0001,0.000137587357554792,0.000189302809589102,0.000260456733490625,0.000358355537183278,0.000493051914261754,0.000678377100206067,0.000933361126430352,0.00128418691029916,0.00176687883594514,0.0024310019015718,0.00334475127847939,0.00460195490083992,0.0063317081439289,0.00871162992331335,0.0119861014114394,0.0164913602058371,0.0226900267320538,0.0312186082091088,0.0429528581002963,0.0590977024544571,0.0813109671826812,0.111873611148895,0.153923945380888,0.211779889093646,0.291382353236417,0.400905280198956,0.551594981323379,0.758924959207687,1.04418479719764,1.4366662704531],\
            [0.0001,0.000137358984481854,0.000188674906178864,0.00025916193509938,0.000355982202216033,0.000488973537900086,0.000671649086044556,0.000922570363872461,0.00126723328294576,0.00174065876847037,0.00239095120766527,0.00328418629830566,0.00451112494784487,0.00619643541706733,0.00851136076295766,0.0116911187095857,0.016058801934055,0.0220582072565705,0.030298929482528,0.0416183018460738,0.0571664767743601,0.0785232919613165,0.107858796419786,0.148153747436569,0.203502482950684,0.279528943976419,0.383958118778863,0.527400972790281,0.724432620372153,0.99507329059848,1.36682256681625],\
            [0.0001,0.000137010983785166,0.00018772009677779,0.000257197151357716,0.000352388347342629,0.000482810741438424,0.00066150374666524,0.000906332791081777,0.00124177547342869,0.00170136879254756,0.00233106212048321,0.0031938111439174,0.00437587206852149,0.0059954253702616,0.00821439128190085,0.0112546183072952,0.0154200632640906,0.0211271803784256,0.0289465576825474,0.0396599634527988,0.054338506095517,0.0744497217756302,0.10200429623011,0.139757089768008,0.191482563600666,0.262352144166329,0.359451253703765,0.492487698927641,0.674762241321687,0.92449838504568,1.26666433242906],\
            [0.0001,0.000136780496792774,0.000187089043028782,0.000255901322499617,0.000350023100214257,0.000478763335362534,0.000654854868570533,0.000895713742502448,0.001225161706836,0.00167578226912512,0.00229214331287458,0.00313520501055222,0.0042883489889053,0.00586562505123261,0.00802303108507742,0.0109739417760076,0.0150102120789731,0.0205310426512685,0.0280824621351416,0.0384113312200896,0.0525392096675567,0.0718633919942818,0.0982951045819178,0.134448532370124,0.183899370506451,0.25153847257751,0.344055572416476,0.470600921194481,0.643690277921186,0.880442759947391,1.20427398103206],\
        ])
    vd_damp_ratio = np.array([\
            [1.68049459715567,1.69724595066663,1.7203522353769,1.75222431023717,1.796187637472,1.85682925399775,1.94047637277048,2.05585655027628,2.21500830058868,2.43453716746209,2.73734831026212,3.15503637793995,3.73118202524743,4.42981740880808,5.12780599874114,6.03610339120981,7.13858359445473,8.51152918174604,9.98408687151498,11.4892724110197,13.0175222001429,14.3659104445271,15.9734565747378,17.206577248012,18.674049267405,19.9324412769751,21.1656704208753,22.1922710080384,23.193307925778,24.1131721623133,24.7553264649226],\
            [1.49805877584555,1.51114965672548,1.52920680181063,1.55411425392304,1.58847080186051,1.6358611323067,1.70122984815926,1.79139737752614,1.91577160220057,2.08732945536551,2.32397090627387,2.65038660368984,3.10063404390184,3.62263200759885,4.28559687496966,5.03693888200916,5.97905744293672,6.76637781624538,7.70207061119717,8.57566387937244,9.78067067253878,10.9288986824501,12.0129474047846,13.3481409969781,14.4109692672984,15.6100488869246,16.6982913511726,17.9267088861263,19.0732946403564,20.0661996643369,21.1718345074578],\
            [0.808635743429545,0.827509858629041,0.853478254993875,0.889207485351604,0.938366389275468,1.00600282618721,1.09906201247829,1.22709968785818,1.40326334208796,1.64564225891467,1.94132779266248,2.27693049703523,2.6361114935881,3.04400433929458,3.5133828328768,4.03083357999218,4.57484357068937,5.17102140871223,5.8984128320865,6.49927334707391,7.18006420150932,7.93498232579605,8.81310434129604,9.73925990088,10.7628298619838,11.9318837828471,13.0316793061949,14.232243788706,15.4549028818545,16.8071925634254,18.3223629569904],\
            [1.08057111007294,1.09127798729422,1.1059848451151,1.12618603566704,1.15393418586244,1.19204876318334,1.24440255953088,1.31631520253156,1.41509367867135,1.55077479038361,1.73714498756528,1.96485910026505,2.26737644370076,2.55306616536831,2.81568923779893,3.17642562310468,3.41309469792577,3.73818091612341,4.18471604413016,4.56620481736329,5.02560433450972,5.56575292637828,6.30769554686197,7.0086699094002,7.91027684925945,8.77450314076682,9.78587142119407,10.7693116847653,11.962670383448,13.2734496122953,14.5764771041722],\
            [1.35081335637868,1.35426387284885,1.35899145941029,1.36546877226742,1.37434340233582,1.38650262029983,1.4031620844529,1.42598738018236,1.45726054241317,1.50010820964633,1.55881422005147,1.63924790244859,1.74758447473341,1.89133035703867,2.06047698218974,2.18861293305607,2.3641732599205,2.60471019089391,2.81892149559473,3.05975220292003,3.3897167242832,3.66106855368764,4.01437862706455,4.43958625986696,4.9898183989673,5.74369686585083,6.4019559727132,7.2844913213172,8.15548287936447,9.24514679065678,10.5036382087092],\
            [1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.14108882860174,1.15406114212932,1.23821599753026,1.3533234268229,1.51076794045477,1.65271047360836,1.72128163014957,1.81507359872319,1.92748795345901,2.06340456071967,2.24931197135468,2.47879134632076,2.68560886127019,2.9684948856725,3.28394065490127,3.66838601511969,4.19423228872321,4.84864489158223,5.59721193650551,6.58821571110552,7.55358756700538,8.63286053423812],\
        ])

    if (PI > np.amax(vd_PI) or (PI < np.amin(vd_PI))):
        print("Invalid PI value for Vucetic-Dobry modulus reduction curve.")
        return -1

    if (PI == np.amax(vd_PI)):
        max_gamma = np.minimum(vd_damp_gamma[-1][-1], vd_MR_gamma[-1][-1])
        gamma = np.logspace(-4, np.log10(max_gamma), 31)
        GGmax_interpolator = interpolate.interp1d(vd_MR_gamma[-1], vd_MR_GGmax[-1])
        damp_interpolator  = interpolate.interp1d(vd_damp_gamma[-1], vd_damp_ratio[-1])
        GGmax = GGmax_interpolator(gamma)
        damp  = damp_interpolator(gamma)
        return GGmax, damp, gamma

    PI_index_interp = interpolate.interp1d(vd_PI, range(len(vd_PI)), kind='previous')
    PI_index = int(PI_index_interp(PI))

    factor_next = (PI - vd_PI[PI_index])/(vd_PI[PI_index+1] - vd_PI[PI_index])
    factor_prev = 1.0 - factor_next

    max_gamma = np.amin(np.array([vd_damp_gamma[PI_index][-1], vd_MR_gamma[PI_index][-1],vd_damp_gamma[PI_index+1][-1], vd_MR_gamma[PI_index+1][-1]]))
    gamma = np.logspace(-4, np.log10(max_gamma), 31)

    GGmax1_interpolator = interpolate.interp1d(vd_MR_gamma[PI_index], vd_MR_GGmax[PI_index])
    damp1_interpolator  = interpolate.interp1d(vd_damp_gamma[PI_index], vd_damp_ratio[PI_index])
    GGmax1 = GGmax1_interpolator(gamma)
    damp1  = damp1_interpolator(gamma)

    GGmax2_interpolator = interpolate.interp1d(vd_MR_gamma[PI_index+1], vd_MR_GGmax[PI_index+1])
    damp2_interpolator  = interpolate.interp1d(vd_damp_gamma[PI_index+1], vd_damp_ratio[PI_index+1])
    GGmax2 = GGmax2_interpolator(gamma)
    damp2  = damp2_interpolator(gamma)

    GGmax = factor_prev * GGmax1 + factor_next * GGmax2
    damp  = factor_prev * damp1 + factor_next * damp2

    # use decimal strains
    gamma /= 100.0

    return GGmax, damp, gamma

def Zhang2005(p, PI, age='Quaternary', Patm=100):
    '''
    This function calculates the Zhange et al. (2005) modulus reduction and 
    damping curves for soil with a confining pressure of 'Conf_pressure', plastic 
    index of 'PI', geologic age of 'AGE', 'Patm' is used to normalize stresses 
    (default Patm = 100 kPa) 

    age = 'Quaternary' | 'Tertiary' | 'Residual'
    
    Outputs are modulus reduction curve 'GGmax', damping ratio (%) 'Damp' and
    corresponding shear strain amplitudes 'gamma'
    
    written by: Alborz Ghofrani
    July 2020 - Golder Associates, Inc.
    '''

    numPoints = 31
    gamma   = np.logspace(-6,-1,numPoints)
    gamma_p = gamma * 100.0

    # constants
    if age == 'Quaternary':
        alpha = 0.0021 * PI + 0.834
        k = 0.310 * np.exp(-0.0142 * PI)
        gamma_r1 = 0.0011 * PI + 0.0749
    elif age == 'Tertiary':
        alpha = 0.0009 * PI + 1.026
        k = 0.316 * np.exp(-0.0110 * PI)
        gamma_r1 = 0.0004 * PI + 0.0311
    elif age == 'Residual':
        alpha = 0.0043 * PI + 0.794
        k = 0.420 * np.exp(-0.0456 * PI)
        gamma_r1 = 0.0009 * PI + 0.0385
    else:
        raise(f'Age {age} is not defined. (Quaternary|Tertiary|Residual)')

    # reference strain (%)
    gr = gamma_r1 * (p / Patm) ** k

    Dmin1 = 0.008 * PI + 0.82
    Dmin = Dmin1 * (p / Patm) ** (-k/2.0)

    GGmax = 1.0 / (1 + (gamma_p/gr)**alpha)
    Damp = Dmin + 10.6 * GGmax**2.0 - 31.6*GGmax + 21.0

    return GGmax, Damp, gamma

def test():
    import matplotlib.pyplot as plt

    plt.figure("vucetic-dobry_MR")
    plt.figure("vucetic-dobry_Damp")
    for PI in np.arange(0,201,15):
        GGmax, damp, gamma = VuceticDobry(PI)
        plt.figure("vucetic-dobry_MR")
        plt.semilogx(gamma, GGmax, label=PI)
        plt.legend()
        plt.grid(which='both')
        plt.figure("vucetic-dobry_Damp")
        plt.semilogx(gamma, damp, label=PI)
        plt.legend()
        plt.grid(which='both')
    plt.show()
    
def test2():
    GGmax, damp, gamma, GGmax_org = Darendeli_Corrected(2400, 0, 4, 2.58e3,2.56,Patm=2117)
    for ii in gamma:
        print(ii)

    for jj in GGmax:
        print(jj)
    
def test3():
    import matplotlib.pyplot as plt

    GGmax, damp, gamma = Zhang2005(200.0, 10.0)

    plt.figure("Zhang_MR")
    plt.semilogx(gamma, GGmax)
    plt.legend()
    plt.grid(which='both')
    plt.figure("Zhang_Damp")
    plt.semilogx(gamma, damp)
    plt.legend()
    plt.grid(which='both')
    plt.show()



if __name__ == "__main__" :
    test3()