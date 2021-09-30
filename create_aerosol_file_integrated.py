import numpy as np
import math
import sys

#READ DATA FROM INPUT FILES
finputs=open("./inputs_params.in",'r')
layer=[]
for i, line in enumerate(finputs):
    if i==13:
        temp = line.split()
        RH = temp[1]
    elif i==11:
        temp = line.split()
        layer.append(temp[1])
    elif i==12:
        temp = line.split()
        layer.append(temp[1])
finputs.close()
if RH=='0':
    RH='00'
fwavelenght=open("./Inputs/integration_limits.dat",'r')
lines=fwavelenght.readlines()
bins=int(lines[0])
for i in range(1,bins+1):
    wl=(float(lines[i])+float(lines[i+1]))/2
    fwavelenght.close()
    for combination_type in layer:
        perc_inso = 0
        perc_waso = 0
        perc_soot = 0
        perc_ssacc = 0
        perc_sscoa = 0
        perc_minnuc = 0
        perc_minacc = 0
        perc_mincoa = 0
        perc_mintrans = 0
        perc_sulf = 0
        perc_fog=0
        if combination_type=='Manual':
            perc_inso = float(input('Particle concentration of insoluble aerosol:'))
            perc_waso = float(input('Particle concentration of water soluble aerosol:'))
            perc_soot = float(input('Particle concentration of soot:'))
            perc_ssacc = float(input('Particle concentration of sea salt (acc):'))
            perc_sscoa = float(input('Particle concentration of sea salt (coa):'))
            perc_minnuc = float(input('Particle concentration of mineral (suc):'))
            perc_minacc = float(input('Particle concentration of mineral (acc):'))
            perc_mincoa = float(input('Particle concentration of mineral (coa):'))
            perc_mintrans = float(input('Particle concentration of mineral transported:'))
            perc_sulf = float(input('Particle concentration of sulfate droplets:'))
            perc_fog = float(input('Particle concentration of fog:'))
        elif combination_type=='CC':
            perc_inso=0.15
            perc_waso=2600
        elif combination_type=='CA':
            perc_inso=0.4
            perc_waso=7000
            perc_soot=8300
        elif combination_type=='CP':
            perc_inso=0.6
            perc_waso=15700
            perc_soot=34300
        elif combination_type == 'U':
            perc_inso=1.5
            perc_waso=28000
            perc_soot=130000
        elif combination_type == 'D':
            perc_waso=2000
            perc_minnuc=269.5
            perc_minacc=30.5
            perc_mincoa=0.142
        elif combination_type == 'MC':
            perc_waso=1500
            perc_ssacc=20
            perc_sscoa=3.2*10**-3
        elif combination_type == 'MP':
            perc_waso=3800
            perc_ssacc=20
            perc_sscoa=3.2*10**-3
            perc_soot=5180
        elif combination_type == 'MT':
            perc_waso=590
            perc_ssacc=10
            perc_sscoa=1.3*10**-3
        elif combination_type == 'ART':
            perc_waso=1300
            perc_inso=0.01
            perc_ssacc=1.9
            perc_soot=5300
        elif combination_type == 'ANT':
            perc_sulf=42.9
            perc_ssacc=0.047
            perc_mintrans=0.0053
        else:
            print('ERROR. Incorret aerosol combination definition. See aerosol_guide.txt')
            sys.exit()

        #obrir els files de cada tipus, buscar wavelength i fer l'average

        ##INSOLUBLE
        finso=open("./Aerosol_optics/OPAC_data/inso00",'r')
        rows, cols = (78-17, 4)
        wl_inso = np.arange(rows*cols, dtype='float')
        wl_inso = wl_inso.reshape(rows, cols)
        pf_inso=np.arange((198-86)*2, dtype='float')
        pf_inso = pf_inso.reshape((198-86), 2)
        a=0
        for i, line in enumerate(finso):
            if i >= 17 and i<78:
                temp = line.split()
                wl_inso[a,0]=float(temp[1]) #wavelength
                wl_inso[a,1] = float(temp[2]) #ext coef
                wl_inso[a,2] = float(temp[3]) #scat coef
                wl_inso[a,3] = float(temp[6]) #ssa param
                a=a+1
            elif i==84:
                temp = line.split()
                temp2 = temp[2:-1]
                absolute_val=[]
                for j in temp2:
                    absolute_val.append(np.abs(float(j) - wl/1000))
                smallest_diff_index = min(range(len(absolute_val)), key=absolute_val.__getitem__)
            if i>=86:
                temp = line.split()
                pf_inso[i-86,0]=temp[0] #scattering angle
                pf_inso[i-86,1]=temp[smallest_diff_index+1]

        finso.close()
        absolute_val_array = np.abs(wl_inso[:,0] - wl/1000)
        smallest_difference_index = absolute_val_array.argmin()
        ext_inso=wl_inso[smallest_difference_index,1]
        scat_inso=wl_inso[smallest_difference_index,2]

        ##SOOT
        fsoot=open("./Aerosol_optics/OPAC_data/soot00",'r')
        rows, cols = (78-17, 4)
        wl_soot = np.arange(rows*cols, dtype='float')
        wl_soot = wl_soot.reshape(rows, cols)
        pf_soot=np.arange((198-86)*2, dtype='float')
        pf_soot = pf_soot.reshape((198-86), 2)
        a=0
        for i, line in enumerate(fsoot):
            if i==17+smallest_difference_index:
                temp = line.split()
                ext_soot=float(temp[2]) #ext coef
                scat_soot=float(temp[3]) #scat coef
            if i >= 86:
                temp = line.split()
                pf_soot[i - 86, 0] = temp[0]  # scattering angle
                pf_soot[i - 86, 1] = temp[smallest_diff_index + 1]
        fsoot.close()

        ##WATER SOLUBLE
        fwaso=open("./Aerosol_optics/OPAC_data/waso"+RH,'r')
        rows, cols = (78-17, 4)
        wl_waso = np.arange(rows*cols, dtype='float')
        wl_waso = wl_waso.reshape(rows, cols)
        pf_waso=np.arange((198-86)*2, dtype='float')
        pf_waso = pf_waso.reshape((198-86), 2)
        for i, line in enumerate(fwaso):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_waso = float(temp[2])  # ext coef
                scat_waso = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_waso[i - 86, 0] = temp[0]  # scattering angle
                pf_waso[i - 86, 1] = temp[smallest_diff_index + 1]
        fwaso.close()

        ##SEA SALT ACC
        fssam=open("./Aerosol_optics/OPAC_data/ssam"+RH,'r')
        rows, cols = (78-17, 4)
        wl_ssam = np.arange(rows*cols, dtype='float')
        wl_ssam = wl_ssam.reshape(rows, cols)
        pf_ssam=np.arange((198-86)*2, dtype='float')
        pf_ssam = pf_ssam.reshape((198-86), 2)
        for i, line in enumerate(fssam):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_ssam = float(temp[2])  # ext coef
                scat_ssam = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_ssam[i - 86, 0] = temp[0]  # scattering angle
                pf_ssam[i - 86, 1] = temp[smallest_diff_index + 1]
        fssam.close()

        ##SEA SALT COA
        fsscm=open("./Aerosol_optics/OPAC_data/sscm"+RH,'r')
        rows, cols = (78-17, 4)
        wl_sscm = np.arange(rows*cols, dtype='float')
        wl_sscm = wl_sscm.reshape(rows, cols)
        pf_sscm=np.arange((198-86)*2, dtype='float')
        pf_sscm = pf_sscm.reshape((198-86), 2)
        for i, line in enumerate(fsscm):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_sscm = float(temp[2])  # ext coef
                scat_sscm = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_sscm[i - 86, 0] = temp[0]  # scattering angle
                pf_sscm[i - 86, 1] = temp[smallest_diff_index + 1]
        fsscm.close()

        ##MINERAL NUC
        fminnuc=open("./Aerosol_optics/OPAC_data/minm00",'r')
        rows, cols = (78-17, 4)
        wl_minnuc = np.arange(rows*cols, dtype='float')
        wl_minnuc = wl_minnuc.reshape(rows, cols)
        pf_minnuc=np.arange((198-86)*2, dtype='float')
        pf_minnuc = pf_minnuc.reshape((198-86), 2)
        for i, line in enumerate(fminnuc):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_minnuc = float(temp[2])  # ext coef
                scat_minnuc = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_minnuc[i - 86, 0] = temp[0]  # scattering angle
                pf_minnuc[i - 86, 1] = temp[smallest_diff_index + 1]
        fminnuc.close()

        ##MINERAL ACC
        fminacc=open("./Aerosol_optics/OPAC_data/miam00",'r')
        rows, cols = (78-17, 4)
        wl_minacc = np.arange(rows*cols, dtype='float')
        wl_minacc = wl_minacc.reshape(rows, cols)
        pf_minacc=np.arange((198-86)*2, dtype='float')
        pf_minacc = pf_minacc.reshape((198-86), 2)
        for i, line in enumerate(fminacc):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_minacc = float(temp[2])  # ext coef
                scat_minacc = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_minacc[i - 86, 0] = temp[0]  # scattering angle
                pf_minacc[i - 86, 1] = temp[smallest_diff_index + 1]
        fminacc.close()

        ##MINERAL COA
        fmincoa=open("./Aerosol_optics/OPAC_data/micm00",'r')
        rows, cols = (78-17, 4)
        wl_mincoa = np.arange(rows*cols, dtype='float')
        wl_mincoa = wl_mincoa.reshape(rows, cols)
        pf_mincoa=np.arange((198-86)*2, dtype='float')
        pf_mincoa = pf_mincoa.reshape((198-86), 2)
        for i, line in enumerate(fmincoa):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_mincoa = float(temp[2])  # ext coef
                scat_mincoa = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_mincoa[i - 86, 0] = temp[0]  # scattering angle
                pf_mincoa[i - 86, 1] = temp[smallest_diff_index + 1]
        fmincoa.close()

        ##MINERAL TRANS
        fmintrans=open("./Aerosol_optics/OPAC_data/mitr00",'r')
        rows, cols = (78-17, 4)
        wl_mintrans = np.arange(rows*cols, dtype='float')
        wl_mintrans = wl_mintrans.reshape(rows, cols)
        pf_mintrans=np.arange((198-86)*2, dtype='float')
        pf_mintrans = pf_mintrans.reshape((198-86), 2)
        for i, line in enumerate(fmintrans):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_mintrans = float(temp[2])  # ext coef
                scat_mintrans = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_mintrans[i - 86, 0] = temp[0]  # scattering angle
                pf_mintrans[i - 86, 1] = temp[smallest_diff_index + 1]
        fmintrans.close()

        ##SULFATE
        fsulf=open("./Aerosol_optics/OPAC_data/suso"+RH,'r')
        rows, cols = (78-17, 4)
        wl_sulf = np.arange(rows*cols, dtype='float')
        wl_sulf = wl_sulf.reshape(rows, cols)
        pf_sulf=np.arange((198-86)*2, dtype='float')
        pf_sulf = pf_sulf.reshape((198-86), 2)
        for i, line in enumerate(fsulf):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_sulf = float(temp[2])  # ext coef
                scat_sulf = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_sulf[i - 86, 0] = temp[0]  # scattering angle
                pf_sulf[i - 86, 1] = temp[smallest_diff_index + 1]
        fsulf.close()

        ##FOG
        ffog=open("./Aerosol_optics/OPAC_data/fogr00",'r')
        rows, cols = (78-17, 4)
        wl_fog = np.arange(rows*cols, dtype='float')
        wl_fog = wl_sulf.reshape(rows, cols)
        pf_fog=np.arange((198-86)*2, dtype='float')
        pf_fog = pf_fog.reshape((198-86), 2)
        for i, line in enumerate(ffog):
            if i == 17 + smallest_difference_index:
                temp = line.split()
                ext_fog = float(temp[2])  # ext coef
                scat_fog = float(temp[3])  # scat coef
            if i >= 86:
                temp = line.split()
                pf_fog[i - 86, 0] = temp[0]  # scattering angle
                pf_fog[i - 86, 1] = temp[smallest_diff_index + 1]
        ffog.close()

        ###crear average de les variables
        #ssa
        ext_total=(ext_inso*perc_inso+ext_soot*perc_soot+ext_waso*perc_waso+perc_ssacc*ext_ssam+perc_sscoa*ext_sscm+perc_minnuc*ext_minnuc+perc_minacc*ext_minacc+perc_mincoa*ext_mincoa+perc_mintrans*ext_mintrans+perc_sulf*ext_sulf+perc_fog*ext_fog)
        scat_total=(scat_inso*perc_inso+scat_soot*perc_soot+scat_waso*perc_waso+perc_ssacc*scat_ssam+perc_sscoa*scat_sscm+perc_minnuc*scat_minnuc+perc_minacc*scat_minacc+perc_mincoa*scat_mincoa+perc_mintrans*scat_mintrans+perc_sulf*scat_sulf+perc_fog*scat_fog)
        w_total=scat_total/ext_total

        #asymmetry
        pf_total=np.arange((198-86)*2, dtype='float')
        pf_total = pf_total.reshape((198-86), 2)
        d=np.arange((198-86)*1, dtype='float')
        d=d.reshape((198-86), 1)

        gtotal_num=0
        gtotal_den=0
        for i in range(len(pf_total)):
            pf_total[i,0]=pf_inso[i,0]*math.pi/180
            pf_total[i,1]=pf_inso[i,1]*perc_inso+pf_soot[i,1]*perc_soot+pf_waso[i,1]*perc_waso+pf_ssam[i,1]*perc_ssacc+pf_sscm[i,1]*perc_sscoa+pf_minnuc[i,1]*perc_minnuc+pf_minacc[i,1]*perc_minacc+pf_mincoa[i,1]*perc_mincoa+pf_mintrans[i,1]*perc_mintrans+pf_sulf[i,1]*perc_sulf+pf_fog[i,1]*perc_fog
        for i in range(len(pf_total)):
            if i==1000:
                gtotal_num = gtotal_num + math.cos(pf_total[i, 0]) * pf_total[i, 1] * (math.cos((pf_total[i+1, 0]-pf_total[i,0])/2) - math.cos(0))
                gtotal_den = gtotal_den + pf_total[i, 1] * (math.cos((pf_total[i+1, 0]-pf_total[i,0])/2) - math.cos(0))
                print(i, pf_total[i, :], gtotal_num, gtotal_den)
            elif i>0  and i<len(pf_total)-1:
                gtotal_num=gtotal_num+math.cos(pf_total[i,0])*pf_total[i,1]*((math.cos(pf_total[i,0]+((pf_total[i+1, 0]-pf_total[i,0])/2)) - math.cos(pf_total[i,0]-((pf_total[i, 0]-pf_total[i-1,0])/2))))
                gtotal_den=gtotal_den+pf_total[i,1]*((math.cos(pf_total[i,0]+((pf_total[i+1, 0]-pf_total[i,0])/2)) - math.cos(pf_total[i,0]-((pf_total[i, 0]-pf_total[i-1,0])/2))))
            elif i==len(pf_total)-1:
                gtotal_num = gtotal_num + math.cos(pf_total[i, 0]) * pf_total[i, 1]*((math.cos(math.pi) - math.cos(pf_total[i,0]-((pf_total[i, 0]-pf_total[i-1,0])/2))))
                gtotal_den = gtotal_den + pf_total[i, 1]*((math.cos(math.pi) - math.cos(pf_total[i,0]-((pf_total[i, 0]-pf_total[i-1,0])/2))))

        g_total=gtotal_num/gtotal_den
        ##normalization

        pf_norm = np.arange((181) * 2, dtype='float')
        pf_norm = pf_norm.reshape((181), 2)
        norm_factor=0

        f=open('./Inputs/'+combination_type+'_'+str(int(wl))+'.txt','w+')
        f.write(str(w_total)+' # single scatering albedo\n')
        f.write('ScatAngle PhaseFct\n')
        for i in range(181):
            pf_norm[i, 0] = i * math.pi / 180
            absolute_val_array = np.abs(pf_total[:, 0]*180/math.pi - i)
            smallest_diffe_ind = absolute_val_array.argmin()
            if min(absolute_val_array)==0:
                smallest_diffe_ind = absolute_val_array.argmin()
                f.write(str(float(i))+' '+str(pf_total[smallest_diffe_ind,1])+'\n')
                pf_norm[i,1]=pf_total[smallest_diffe_ind,1]
            else:
                if pf_total[smallest_diffe_ind,0]*180/math.pi - i < 0:
                    pf=(pf_total[smallest_diffe_ind,1]*np.abs(1/(pf_total[smallest_diffe_ind,0]*180/math.pi - i))+pf_total[smallest_diffe_ind+1,1]*np.abs(1/(pf_total[smallest_diffe_ind+1,0]*180/math.pi - i)))/(np.abs(1/(pf_total[smallest_diffe_ind,0]*180/math.pi - i))+np.abs(1/(pf_total[smallest_diffe_ind+1,0]*180/math.pi - i)))
                else:
                    pf = (pf_total[smallest_diffe_ind, 1] * np.abs(1 / (pf_total[smallest_diffe_ind, 0] * 180 / math.pi - i)) +
                          pf_total[smallest_diffe_ind - 1, 1] * np.abs(1 / (pf_total[smallest_diffe_ind - 1, 0] * 180 / math.pi - i))) / (
                                     np.abs(1 / (pf_total[smallest_diffe_ind, 0] * 180 / math.pi - i)) + np.abs(
                                 1 / (pf_total[smallest_diffe_ind - 1, 0] * 180 / math.pi - i)))
                f.write(str(float(i)) + ' ' + str(pf) + '\n')
                pf_norm[i, 1] = pf
        #print(w_total,g_total)
        f.close()

        pf_norm2 = np.arange((181) * 2, dtype='float')
        pf_norm2 = pf_norm2.reshape((181), 2)
        norm = 1
        if norm == 1:
            for i in range(181):
                norm_factor = norm_factor+ (1 / (4 * math.pi)) * 2 * math.pi * pf_norm[i, 1] * math.sin(pf_norm[i, 0]) * math.pi / 180
            for i in range(181):
                pf_norm2[i, 1] = pf_norm[i, 1] / norm_factor
                pf_norm2[i, 0] = pf_norm[i, 0]

        print(pf_norm)
        norm_factor2=0
        g_total_norm=0
        for i in range(181):
            norm_factor2 =norm_factor2 + (1 / (4 * math.pi)) * 2 * math.pi * pf_norm[i, 1] * math.sin(pf_norm[i, 0]) * math.pi / 180
            g_total_norm=g_total_norm+(0.5*math.sin(pf_norm[i,0])*pf_norm2[i,1]*math.cos(pf_norm[i,0])*math.pi/180)
        print(norm_factor2,g_total_norm)
        

            
