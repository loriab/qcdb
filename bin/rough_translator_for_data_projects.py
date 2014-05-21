
#with open('../../fullcurve/gath_igor_reapwaves.dat', 'r') as f:
with open('accumulate', 'r') as f:
    content = f.readlines()

f1 = open('S22_dhdft.py', 'w')
f2 = open('NBC1_dhdft.py', 'w')
f3 = open('HBC1_dhdft.py', 'w')
f4 = open('HSG_dhdft.py', 'w')
f1.write('\ndef load_dhdft(dbinstance):\n\n')
f2.write('\ndef load_dhdft(dbinstance):\n\n')
f3.write('\ndef load_dhdft(dbinstance):\n\n')
f4.write('\ndef load_dhdft(dbinstance):\n\n')


for line in content:

    # skip separations btwn modelchems
    if line == '\n':
        continue

    # read modelchem for following entries
    lline = line.strip().split()
    if lline[0].endswith('-rcurve'):
        #database, method, basis, tmp1, tmp2 = lline[2].split('-')
        mtdarr = lline[2].split('-')
        if mtdarr[2] == 'CP':
            database = mtdarr[0]
            method = mtdarr[1]
            bsse = mtdarr[2]
            basis = mtdarr[3]
            role = mtdarr[4]
        else:
            database = mtdarr[0]
            method = mtdarr[1]
            bsse = 'unCP'
            basis = mtdarr[2]
            role = mtdarr[3]

#S22/S22-B3LYPD2-CP-adz/reap_S22-B3LYPD2-CP-adz-DASH.outie
    #S22/S22-B3LYPD2-CP-adz/reap_S22-B3LYPD2-CP-adz-RAW.outie
#S22/S22-B3LYPD2-CP-atz/reap_S22-B3LYPD2-CP-atz-DASH.outie
    #S22/S22-B3LYPD2-CP-atz/reap_S22-B3LYPD2-CP-atz-RAW.outie
#S22/S22-dlDFD-adz/reap_S22-dlDFD-adz-DASH.outie
    #S22/S22-dlDFD-adz/reap_S22-dlDFD-adz-RAW.outie
#S22/S22-dlDFD-atz/reap_S22-dlDFD-atz-DASH.outie
    #S22/S22-dlDFD-atz/reap_S22-dlDFD-atz-RAW.outie
#S22/S22-dlDFD-CP-adz/reap_S22-dlDFD-CP-adz-DASH.outie
#S22/S22-dlDFD-CP-adz/reap_S22-dlDFD-CP-adz-RAW.outie
#S22/S22-dlDFD-CP-atz/reap_S22-dlDFD-CP-atz-DASH.outie
#S22/S22-dlDFD-CP-atz/reap_S22-dlDFD-CP-atz-RAW.outie
#S22/S22-DSDPBEP86-adz/reap_S22-DSDPBEP86-adz-RAW.outie
#S22/S22-DSDPBEP86-atz/reap_S22-DSDPBEP86-atz-RAW.outie
#S22/S22-DSDPBEP86-CP-adz/reap_S22-DSDPBEP86-CP-adz-RAW.outie
#S22/S22-DSDPBEP86-CP-atz/reap_S22-DSDPBEP86-CP-atz-RAW.outie
#S22/S22-LCVV10-adz/reap_S22-LCVV10-adz-RAW.outie
#S22/S22-LCVV10-atz/reap_S22-LCVV10-atz-RAW.outie
#S22/S22-LCVV10-CP-adz/reap_S22-LCVV10-CP-adz-RAW.outie
#S22/S22-LCVV10-CP-atz/reap_S22-LCVV10-CP-atz-RAW.outie
#S22/S22-PBE02-adz/reap_S22-PBE02-adz-RAW.outie
#S22/S22-PBE02-atz/reap_S22-PBE02-atz-RAW.outie
#S22/S22-PBE02-CP-adz/reap_S22-PBE02-CP-adz-RAW.outie
#S22/S22-PBE02-CP-atz/reap_S22-PBE02-CP-atz-RAW.outie
#S22/S22-VV10-adz/reap_S22-VV10-adz-RAW.outie
#S22/S22-VV10-atz/reap_S22-VV10-atz-RAW.outie
#S22/S22-VV10-CP-adz/reap_S22-VV10-CP-adz-RAW.outie
#S22/S22-VV10-CP-atz/reap_S22-VV10-CP-atz-RAW.outie
#S22/S22-wB97X2-adz/reap_S22-wB97X2-adz-RAW.outie
#S22/S22-wB97X2-atz/reap_S22-wB97X2-atz-RAW.outie
#S22/S22-wB97X2-CP-adz/reap_S22-wB97X2-CP-adz-RAW.outie
#S22/S22-wB97X2-CP-atz/reap_S22-wB97X2-CP-atz-RAW.outie

    # read all the good data
    else:
        if not ((method == 'B3LYPD2' and role == 'RAW') or (method == 'dlDFD' and role == 'RAW')):
            if lline[0] == '1':
                print database, method, basis, bsse, role
            if database == 'S22':
                f1.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" % 
                    (database, lline[0], method, bsse, basis, float(lline[2])))
            elif database == 'NBC1':
                f2.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" % 
                    (database, lline[0], method, bsse, basis, float(lline[2])))
            elif database == 'HBC1':
                f3.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" % 
                    (database, lline[0], method, bsse, basis, float(lline[2])))
            elif database == 'HSG':
                f4.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" % 
                    (database, lline[0], method, bsse, basis, float(lline[2])))

f1.write('\n')
f2.write('\n')
f3.write('\n')
f4.write('\n')
f1.close()
f2.close()
f3.close()
f4.close()
