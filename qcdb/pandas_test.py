

def test_pandas(h2kc, project, df):
    """

    """
    import qcdb

    try:

        if project == 'dhdft':
            print """\n        Checks from dhdft project:\n"""
            digits = 4

            # B3LYP
            qcdb.compare_values(-2.8732, h2kc * df['B3LYP-CP-adz']['S22-22'], digits, 'B3LYP-CP-adz')
            qcdb.compare_values(-2.9389, h2kc * df['B3LYP-CP-atz']['S22-22'], digits, 'B3LYP-CP-atz')
            qcdb.compare_values(-3.6276, h2kc * df['B3LYP-unCP-adz']['S22-22'], digits, 'B3LYP-unCP-adz')
            qcdb.compare_values(-3.1766, h2kc * df['B3LYP-unCP-atz']['S22-22'], digits, 'B3LYP-unCP-atz')
            qcdb.compare_values(-7.1256, h2kc * df['B3LYPD3-CP-adz']['S22-22'], digits, 'B3LYPD3-CP-adz')
            qcdb.compare_values(-7.1913, h2kc * df['B3LYPD3-CP-atz']['S22-22'], digits, 'B3LYPD3-CP-atz')
            qcdb.compare_values(-7.8800, h2kc * df['B3LYPD3-unCP-adz']['S22-22'], digits, 'B3LYPD3-unCP-adz')
            qcdb.compare_values(-7.4290, h2kc * df['B3LYPD3-unCP-atz']['S22-22'], digits, 'B3LYPD3-unCP-atz')

            # B2PLYP
            qcdb.compare_values(-4.5423, h2kc * df['B2PLYP-nfc-CP-adz']['S22-22'], digits, 'B2PLYP-nfc-CP-adz')
            qcdb.compare_values(-4.7537, h2kc * df['B2PLYP-nfc-CP-atz']['S22-22'], digits, 'B2PLYP-nfc-CP-atz')
            qcdb.compare_values(-6.1559, h2kc * df['B2PLYP-nfc-unCP-adz']['S22-22'], digits, 'B2PLYP-nfc-unCP-adz')
            qcdb.compare_values(-5.8908, h2kc * df['B2PLYP-nfc-unCP-atz']['S22-22'], digits, 'B2PLYP-nfc-unCP-atz')
            qcdb.compare_values(-6.7902, h2kc * df['B2PLYPD3-nfc-CP-adz']['S22-22'], digits, 'B2PLYPD3-nfc-CP-adz')
            qcdb.compare_values(-7.0016, h2kc * df['B2PLYPD3-nfc-CP-atz']['S22-22'], digits, 'B2PLYPD3-nfc-CP-atz')
            qcdb.compare_values(-8.4038, h2kc * df['B2PLYPD3-nfc-unCP-adz']['S22-22'], digits, 'B2PLYPD3-nfc-unCP-adz')
            qcdb.compare_values(-8.1387, h2kc * df['B2PLYPD3-nfc-unCP-atz']['S22-22'], digits, 'B2PLYPD3-nfc-unCP-atz')

            # B97D
            qcdb.compare_values(-6.2894, h2kc * df['B97D3-CP-adz']['S22-22'], digits, 'B97D3-CP-adz')
            qcdb.compare_values(-6.3362, h2kc * df['B97D3-CP-atz']['S22-22'], digits, 'B97D3-CP-atz')
            qcdb.compare_values(-7.0221, h2kc * df['B97D3-unCP-adz']['S22-22'], digits, 'B97D3-unCP-adz')
            qcdb.compare_values(-6.5884, h2kc * df['B97D3-unCP-atz']['S22-22'], digits, 'B97D3-unCP-atz')

            # M052X & M062X
            qcdb.compare_values(-5.9540, h2kc * df['M052X-CP-adz']['S22-22'], digits, 'M052X-CP-adz')
            qcdb.compare_values(-6.0720, h2kc * df['M052X-CP-atz']['S22-22'], digits, 'M052X-CP-atz')
            qcdb.compare_values(-6.7825, h2kc * df['M052X-unCP-adz']['S22-22'], digits, 'M052X-unCP-adz')
            qcdb.compare_values(-6.3942, h2kc * df['M052X-unCP-atz']['S22-22'], digits, 'M052X-unCP-atz')
            qcdb.compare_values(-6.5146, h2kc * df['M062X-CP-adz']['S22-22'], digits, 'M062X-CP-adz')
            qcdb.compare_values(-6.5581, h2kc * df['M062X-CP-atz']['S22-22'], digits, 'M062X-CP-atz')
            qcdb.compare_values(-7.3393, h2kc * df['M062X-unCP-adz']['S22-22'], digits, 'M062X-unCP-adz')
            qcdb.compare_values(-6.8876, h2kc * df['M062X-unCP-atz']['S22-22'], digits, 'M062X-unCP-atz')

            # VV10 & LCVV10
            qcdb.compare_values(-6.9325, h2kc * df['VV10-CP-adz']['S22-22'], digits, 'VV10-CP-adz')
            qcdb.compare_values(-6.9913, h2kc * df['VV10-CP-atz']['S22-22'], digits, 'VV10-CP-atz')
            qcdb.compare_values(-7.7285, h2kc * df['VV10-unCP-adz']['S22-22'], digits, 'VV10-unCP-adz')
            qcdb.compare_values(-7.2664, h2kc * df['VV10-unCP-atz']['S22-22'], digits, 'VV10-unCP-atz')
            qcdb.compare_values(-6.9236, h2kc * df['LCVV10-CP-adz']['S22-22'], digits, 'LCVV10-CP-adz')
            qcdb.compare_values(-6.9086, h2kc * df['LCVV10-CP-atz']['S22-22'], digits, 'LCVV10-CP-atz')
            qcdb.compare_values(-7.7457, h2kc * df['LCVV10-unCP-adz']['S22-22'], digits, 'LCVV10-unCP-adz')
            qcdb.compare_values(-7.1775, h2kc * df['LCVV10-unCP-atz']['S22-22'], digits, 'LCVV10-unCP-atz')

            # DSD-PBEP86
            qcdb.compare_values(-5.4858, h2kc * df['DSDPBEP86-nfc-CP-adz']['S22-22'], digits, 'DSDPBEP86-nfc-CP-adz')
            qcdb.compare_values(-5.7625, h2kc * df['DSDPBEP86-nfc-CP-atz']['S22-22'], digits, 'DSDPBEP86-nfc-CP-atz')
            qcdb.compare_values(-7.6463, h2kc * df['DSDPBEP86-nfc-unCP-adz']['S22-22'], digits, 'DSDPBEP86-nfc-unCP-adz')
            qcdb.compare_values(-7.5912, h2kc * df['DSDPBEP86-nfc-unCP-atz']['S22-22'], digits, 'DSDPBEP86-nfc-unCP-atz')
            # TODO need -D !!

            # PBE & PBE0 & PBE0-2
            qcdb.compare_values(-3.8210, h2kc * df['PBE-CP-adz']['S22-22'], digits, 'PBE-CP-adz')
            qcdb.compare_values(-3.8560, h2kc * df['PBE-CP-atz']['S22-22'], digits, 'PBE-CP-atz')
            qcdb.compare_values(-4.5668, h2kc * df['PBE-unCP-adz']['S22-22'], digits, 'PBE-unCP-adz')
            qcdb.compare_values(-4.0963, h2kc * df['PBE-unCP-atz']['S22-22'], digits, 'PBE-unCP-atz')
            qcdb.compare_values(-4.1464, h2kc * df['PBE0-CP-adz']['S22-22'], digits, 'PBE0-CP-adz')
            qcdb.compare_values(-4.1539, h2kc * df['PBE0-CP-atz']['S22-22'], digits, 'PBE0-CP-atz')
            qcdb.compare_values(-4.9312, h2kc * df['PBE0-unCP-adz']['S22-22'], digits, 'PBE0-unCP-adz')
            qcdb.compare_values(-4.4044, h2kc * df['PBE0-unCP-atz']['S22-22'], digits, 'PBE0-unCP-atz')
            qcdb.compare_values(-6.9963, h2kc * df['PBE0D3-CP-adz']['S22-22'], digits, 'PBE0D3-CP-adz')
            qcdb.compare_values(-7.0039, h2kc * df['PBE0D3-CP-atz']['S22-22'], digits, 'PBE0D3-CP-atz')
            qcdb.compare_values(-7.7811, h2kc * df['PBE0D3-unCP-adz']['S22-22'], digits, 'PBE0D3-unCP-adz')
            qcdb.compare_values(-7.2544, h2kc * df['PBE0D3-unCP-atz']['S22-22'], digits, 'PBE0D3-unCP-atz')
            qcdb.compare_values(-6.6748, h2kc * df['PBED3-CP-adz']['S22-22'], digits, 'PBED3-CP-adz')
            qcdb.compare_values(-6.7098, h2kc * df['PBED3-CP-atz']['S22-22'], digits, 'PBED3-CP-atz')
            qcdb.compare_values(-7.4206, h2kc * df['PBED3-unCP-adz']['S22-22'], digits, 'PBED3-unCP-adz')
            qcdb.compare_values(-6.9501, h2kc * df['PBED3-unCP-atz']['S22-22'], digits, 'PBED3-unCP-atz')
            qcdb.compare_values(-6.0202, h2kc * df['PBE02-CP-adz']['S22-22'], digits, 'PBE02-CP-adz')
            qcdb.compare_values(-6.2983, h2kc * df['PBE02-CP-atz']['S22-22'], digits, 'PBE02-CP-atz')
            qcdb.compare_values(-8.3761, h2kc * df['PBE02-unCP-adz']['S22-22'], digits, 'PBE02-unCP-adz')
            qcdb.compare_values(-8.1894, h2kc * df['PBE02-unCP-atz']['S22-22'], digits, 'PBE02-unCP-atz')

            # higher M0N
            qcdb.compare_values(-6.5554, h2kc * df['M08HX-CP-adz']['S22-22'], digits, 'M08HX-CP-adz')
            qcdb.compare_values(-6.3784, h2kc * df['M08HX-CP-atz']['S22-22'], digits, 'M08HX-CP-atz')
            qcdb.compare_values(-7.5260, h2kc * df['M08HX-unCP-adz']['S22-22'], digits, 'M08HX-unCP-adz')
            qcdb.compare_values(-6.9478, h2kc * df['M08HX-unCP-atz']['S22-22'], digits, 'M08HX-unCP-atz')
            qcdb.compare_values(-6.6251, h2kc * df['M08SO-CP-adz']['S22-22'], digits, 'M08SO-CP-adz')
            qcdb.compare_values(-6.2753, h2kc * df['M08SO-CP-atz']['S22-22'], digits, 'M08SO-CP-atz')
            qcdb.compare_values(-7.5524, h2kc * df['M08SO-unCP-adz']['S22-22'], digits, 'M08SO-unCP-adz')
            qcdb.compare_values(-6.8447, h2kc * df['M08SO-unCP-atz']['S22-22'], digits, 'M08SO-unCP-atz')
            qcdb.compare_values(-6.4102, h2kc * df['M11-CP-adz']['S22-22'], digits, 'M11-CP-adz')
            qcdb.compare_values(-6.0122, h2kc * df['M11-CP-atz']['S22-22'], digits, 'M11-CP-atz')
            qcdb.compare_values(-7.4823, h2kc * df['M11-unCP-adz']['S22-22'], digits, 'M11-unCP-adz')
            qcdb.compare_values(-6.8535, h2kc * df['M11-unCP-atz']['S22-22'], digits, 'M11-unCP-atz')
            qcdb.compare_values(-4.4887, h2kc * df['M11L-CP-adz']['S22-22'], digits, 'M11L-CP-adz')
            qcdb.compare_values(-4.6106, h2kc * df['M11L-CP-atz']['S22-22'], digits, 'M11L-CP-atz')
            qcdb.compare_values(-5.8837, h2kc * df['M11L-unCP-adz']['S22-22'], digits, 'M11L-unCP-adz')
            qcdb.compare_values(-5.4856, h2kc * df['M11L-unCP-atz']['S22-22'], digits, 'M11L-unCP-atz')

            # WB97X-D & WB97X-2
            qcdb.compare_values(-6.9619, h2kc * df['WB97XD-CP-adz']['S22-22'], digits, 'WB97XD-CP-adz')
            qcdb.compare_values(-6.9079, h2kc * df['WB97XD-CP-atz']['S22-22'], digits, 'WB97XD-CP-atz')
            qcdb.compare_values(-7.7526, h2kc * df['WB97XD-unCP-adz']['S22-22'], digits, 'WB97XD-unCP-adz')
            qcdb.compare_values(-7.2837, h2kc * df['WB97XD-unCP-atz']['S22-22'], digits, 'WB97XD-unCP-atz')
            qcdb.compare_values(-6.6073, h2kc * df['WB97X2-nfc-CP-adz']['S22-22'], digits, 'WB97X2-nfc-CP-adz')
            qcdb.compare_values(-6.9328, h2kc * df['WB97X2-nfc-CP-atz']['S22-22'], digits, 'WB97X2-nfc-CP-atz')
            qcdb.compare_values(-8.7650, h2kc * df['WB97X2-nfc-unCP-adz']['S22-22'], digits, 'WB97X2-nfc-unCP-adz')
            qcdb.compare_values(-8.6417, h2kc * df['WB97X2-nfc-unCP-atz']['S22-22'], digits, 'WB97X2-nfc-unCP-atz')

            # dlDFD
            qcdb.compare_values(-6.3829, h2kc * df['DLDFD-CP-adz']['S22-22'], digits, 'DLDFD-CP-adz')
            qcdb.compare_values(-6.4839, h2kc * df['DLDFD-CP-atz']['S22-22'], digits, 'DLDFD-CP-atz')
            qcdb.compare_values(-7.2422, h2kc * df['DLDFD-unCP-adz']['S22-22'], digits, 'DLDFD-unCP-adz')
            qcdb.compare_values(-6.8258, h2kc * df['DLDFD-unCP-atz']['S22-22'], digits, 'DLDFD-unCP-atz')

            # duplicates of checks from dft project that should match dhdft project
            #   digits shouldn't need to be this big but is: df vs non-df & psi4 vs qchem
            print """\n        Checks from dft project:\n"""
            digits = 0.04

            # B3LYP
            qcdb.compare_values(-3.6463, h2kc * df['B3LYP-unCP-adz']['S22-22'], digits, 'B3LYP-unCP-adz')
            qcdb.compare_values(-3.1681, h2kc * df['B3LYP-unCP-atz']['S22-22'], digits, 'B3LYP-unCP-atz')
            qcdb.compare_values(-7.8987, h2kc * df['B3LYPD3-unCP-adz']['S22-22'], digits, 'B3LYPD3-unCP-adz')
            qcdb.compare_values(-7.4205, h2kc * df['B3LYPD3-unCP-atz']['S22-22'], digits, 'B3LYPD3-unCP-atz')
            qcdb.compare_values(-2.8746, h2kc * df['B3LYP-CP-adz']['S22-22'], digits, 'B3LYP-CP-adz')
            qcdb.compare_values(-2.9407, h2kc * df['B3LYP-CP-atz']['S22-22'], digits, 'B3LYP-CP-atz')
            qcdb.compare_values(-7.1270, h2kc * df['B3LYPD3-CP-adz']['S22-22'], digits, 'B3LYPD3-CP-adz')
            qcdb.compare_values(-7.1931, h2kc * df['B3LYPD3-CP-atz']['S22-22'], digits, 'B3LYPD3-CP-atz')

            # B2PLYP
            # none match b/c core not frozen in dhdft project

            # B97-D
            qcdb.compare_values(-7.0599, h2kc * df['B97D3-unCP-adz']['S22-22'], digits, 'B97D3-unCP-adz')
            qcdb.compare_values(-6.5760, h2kc * df['B97D3-unCP-atz']['S22-22'], digits, 'B97D3-unCP-atz')

            # wB97X-D
            qcdb.compare_values(-6.9661, h2kc * df['WB97XD-CP-adz']['S22-22'], digits, 'WB97XD-CP-adz')
            qcdb.compare_values(-6.9069, h2kc * df['WB97XD-CP-atz']['S22-22'], digits, 'WB97XD-CP-atz')
            qcdb.compare_values(-7.7283, h2kc * df['WB97XD-unCP-adz']['S22-22'], digits, 'WB97XD-unCP-adz')
            #qcdb.compare_values(-7.2139, h2kc * df['WB97XD-unCP-atz']['S22-22'], digits, 'WB97XD-unCP-atz')
            # aggh, CP are right on but unCP-adz off by 0.02 and unCP-atz by 0.06!

            # M05-2X
            qcdb.compare_values(-5.9556, h2kc * df['M052X-CP-adz']['S22-22'], digits, 'M052X-CP-adz')
            qcdb.compare_values(-6.0724, h2kc * df['M052X-CP-atz']['S22-22'], digits, 'M052X-CP-atz')
            qcdb.compare_values(-6.7954, h2kc * df['M052X-unCP-adz']['S22-22'], digits, 'M052X-unCP-adz')
            qcdb.compare_values(-6.4217, h2kc * df['M052X-unCP-atz']['S22-22'], digits, 'M052X-unCP-atz')

            # M06-2X
            qcdb.compare_values(-6.5146, h2kc * df['M062X-CP-adz']['S22-22'], digits, 'M062X-CP-adz')
            qcdb.compare_values(-6.5581, h2kc * df['M062X-CP-atz']['S22-22'], digits, 'M062X-CP-atz')
            qcdb.compare_values(-7.3393, h2kc * df['M062X-unCP-adz']['S22-22'], digits, 'M062X-unCP-adz')
            qcdb.compare_values(-6.8876, h2kc * df['M062X-unCP-atz']['S22-22'], digits, 'M062X-unCP-atz')
            qcdb.compare_values(-7.8457, h2kc * df['M062XD3-unCP-adz']['S22-22'], digits, 'M062XD3-unCP-adz')
            qcdb.compare_values(-7.3941, h2kc * df['M062XD3-unCP-atz']['S22-22'], digits, 'M062XD3-unCP-atz')

            # PBE
            qcdb.compare_values(-3.8225, h2kc * df['PBE-CP-adz']['S22-22'], digits, 'PBE-CP-adz')
            qcdb.compare_values(-3.8524, h2kc * df['PBE-CP-atz']['S22-22'], digits, 'PBE-CP-atz')
            qcdb.compare_values(-4.5854, h2kc * df['PBE-unCP-adz']['S22-22'], digits, 'PBE-unCP-adz')
            qcdb.compare_values(-4.0802, h2kc * df['PBE-unCP-atz']['S22-22'], digits, 'PBE-unCP-atz')
            qcdb.compare_values(-6.6763, h2kc * df['PBED3-CP-adz']['S22-22'], digits, 'PBED3-CP-adz')
            qcdb.compare_values(-6.7063, h2kc * df['PBED3-CP-atz']['S22-22'], digits, 'PBED3-CP-atz')
            qcdb.compare_values(-7.4392, h2kc * df['PBED3-unCP-adz']['S22-22'], digits, 'PBED3-unCP-adz')
            qcdb.compare_values(-6.9340, h2kc * df['PBED3-unCP-atz']['S22-22'], digits, 'PBED3-unCP-atz')

            # PBE0
            qcdb.compare_values(-4.9388, h2kc * df['PBE0-unCP-adz']['S22-22'], digits, 'PBE0-unCP-adz')
            qcdb.compare_values(-4.3888, h2kc * df['PBE0-unCP-atz']['S22-22'], digits, 'PBE0-unCP-atz')
            qcdb.compare_values(-7.7887, h2kc * df['PBE0D3-unCP-adz']['S22-22'], digits, 'PBE0D3-unCP-adz')  # new dash
            qcdb.compare_values(-7.2387, h2kc * df['PBE0D3-unCP-atz']['S22-22'], digits, 'PBE0D3-unCP-atz')  # new dash

        elif project == 'dft':
            print """\n        Checks from dft project:\n"""
            digits = 4

            # B3LYP
            qcdb.compare_values(-3.6463, h2kc * df['B3LYP-unCP-adz']['S22-22'], digits, 'B3LYP-unCP-adz')
            qcdb.compare_values(-3.1681, h2kc * df['B3LYP-unCP-atz']['S22-22'], digits, 'B3LYP-unCP-atz')
            qcdb.compare_values(-7.9353, h2kc * df['B3LYPD2-unCP-adz']['S22-22'], digits, 'B3LYPD2-unCP-adz')
            qcdb.compare_values(-7.8987, h2kc * df['B3LYPD3-unCP-adz']['S22-22'], digits, 'B3LYPD3-unCP-adz')
            qcdb.compare_values(-7.4571, h2kc * df['B3LYPD2-unCP-atz']['S22-22'], digits, 'B3LYPD2-unCP-atz')
            qcdb.compare_values(-7.4205, h2kc * df['B3LYPD3-unCP-atz']['S22-22'], digits, 'B3LYPD3-unCP-atz')
            qcdb.compare_values(-2.8746, h2kc * df['B3LYP-CP-adz']['S22-22'], digits, 'B3LYP-CP-adz')
            qcdb.compare_values(-2.9407, h2kc * df['B3LYP-CP-atz']['S22-22'], digits, 'B3LYP-CP-atz')
            qcdb.compare_values(-7.1636, h2kc * df['B3LYPD2-CP-adz']['S22-22'], digits, 'B3LYPD2-CP-adz')
            qcdb.compare_values(-7.1270, h2kc * df['B3LYPD3-CP-adz']['S22-22'], digits, 'B3LYPD3-CP-adz')
            qcdb.compare_values(-7.2297, h2kc * df['B3LYPD2-CP-atz']['S22-22'], digits, 'B3LYPD2-CP-atz')
            qcdb.compare_values(-7.1931, h2kc * df['B3LYPD3-CP-atz']['S22-22'], digits, 'B3LYPD3-CP-atz')
            qcdb.compare_values(-2.8917, h2kc * df['B3LYP-CP-dadz']['S22-22'], digits, 'B3LYP-CP-dadz')
            qcdb.compare_values(-2.9393, h2kc * df['B3LYP-CP-datz']['S22-22'], digits, 'B3LYP-CP-datz')
            qcdb.compare_values(-3.0227, h2kc * df['B3LYP-CP-dz']['S22-22'], digits, 'B3LYP-CP-dz')
            qcdb.compare_values(-2.8620, h2kc * df['B3LYP-CP-hadz']['S22-22'], digits, 'B3LYP-CP-hadz')
            qcdb.compare_values(-2.9278, h2kc * df['B3LYP-CP-hatz']['S22-22'], digits, 'B3LYP-CP-hatz')
            qcdb.compare_values(-2.9521, h2kc * df['B3LYP-CP-tz']['S22-22'], digits, 'B3LYP-CP-tz')
            qcdb.compare_values(-3.6463, h2kc * df['B3LYP-unCP-adz']['S22-22'], digits, 'B3LYP-unCP-adz')
            qcdb.compare_values(-3.1681, h2kc * df['B3LYP-unCP-atz']['S22-22'], digits, 'B3LYP-unCP-atz')
            qcdb.compare_values(-3.9095, h2kc * df['B3LYP-unCP-dadz']['S22-22'], digits, 'B3LYP-unCP-dadz')
            qcdb.compare_values(-3.1769, h2kc * df['B3LYP-unCP-datz']['S22-22'], digits, 'B3LYP-unCP-datz')
            qcdb.compare_values(-6.3591, h2kc * df['B3LYP-unCP-dz']['S22-22'], digits, 'B3LYP-unCP-dz')
            qcdb.compare_values(-3.5103, h2kc * df['B3LYP-unCP-hadz']['S22-22'], digits, 'B3LYP-unCP-hadz')
            qcdb.compare_values(-3.1173, h2kc * df['B3LYP-unCP-hatz']['S22-22'], digits, 'B3LYP-unCP-hatz')
            qcdb.compare_values(-3.3748, h2kc * df['B3LYP-unCP-6311pg_3df_2p_']['S22-22'], digits, 'B3LYP-unCP-6311pg_3df_2p_')
            qcdb.compare_values(-4.1382, h2kc * df['B3LYP-unCP-tz']['S22-22'], digits, 'B3LYP-unCP-tz')
            qcdb.compare_values(-7.1636, h2kc * df['B3LYPD2-CP-adz']['S22-22'], digits, 'B3LYPD2-CP-adz')
            qcdb.compare_values(-7.1270, h2kc * df['B3LYPD3-CP-adz']['S22-22'], digits, 'B3LYPD3-CP-adz')
            qcdb.compare_values(-7.2297, h2kc * df['B3LYPD2-CP-atz']['S22-22'], digits, 'B3LYPD2-CP-atz')
            qcdb.compare_values(-7.1931, h2kc * df['B3LYPD3-CP-atz']['S22-22'], digits, 'B3LYPD3-CP-atz')
            qcdb.compare_values(-7.1807, h2kc * df['B3LYPD2-CP-dadz']['S22-22'], digits, 'B3LYPD2-CP-dadz')
            qcdb.compare_values(-7.1441, h2kc * df['B3LYPD3-CP-dadz']['S22-22'], digits, 'B3LYPD3-CP-dadz')
            qcdb.compare_values(-7.2284, h2kc * df['B3LYPD2-CP-datz']['S22-22'], digits, 'B3LYPD2-CP-datz')
            qcdb.compare_values(-7.1917, h2kc * df['B3LYPD3-CP-datz']['S22-22'], digits, 'B3LYPD3-CP-datz')
            qcdb.compare_values(-7.3117, h2kc * df['B3LYPD2-CP-dz']['S22-22'], digits, 'B3LYPD2-CP-dz')
            qcdb.compare_values(-7.2751, h2kc * df['B3LYPD3-CP-dz']['S22-22'], digits, 'B3LYPD3-CP-dz')
            qcdb.compare_values(-7.1510, h2kc * df['B3LYPD2-CP-hadz']['S22-22'], digits, 'B3LYPD2-CP-hadz')
            qcdb.compare_values(-7.1144, h2kc * df['B3LYPD3-CP-hadz']['S22-22'], digits, 'B3LYPD3-CP-hadz')
            qcdb.compare_values(-7.2168, h2kc * df['B3LYPD2-CP-hatz']['S22-22'], digits, 'B3LYPD2-CP-hatz')
            qcdb.compare_values(-7.1802, h2kc * df['B3LYPD3-CP-hatz']['S22-22'], digits, 'B3LYPD3-CP-hatz')
            qcdb.compare_values(-7.2411, h2kc * df['B3LYPD2-CP-tz']['S22-22'], digits, 'B3LYPD2-CP-tz')
            qcdb.compare_values(-7.2045, h2kc * df['B3LYPD3-CP-tz']['S22-22'], digits, 'B3LYPD3-CP-tz')
            qcdb.compare_values(-7.9353, h2kc * df['B3LYPD2-unCP-adz']['S22-22'], digits, 'B3LYPD2-unCP-adz')
            qcdb.compare_values(-7.8987, h2kc * df['B3LYPD3-unCP-adz']['S22-22'], digits, 'B3LYPD3-unCP-adz')
            qcdb.compare_values(-7.4571, h2kc * df['B3LYPD2-unCP-atz']['S22-22'], digits, 'B3LYPD2-unCP-atz')
            qcdb.compare_values(-7.4205, h2kc * df['B3LYPD3-unCP-atz']['S22-22'], digits, 'B3LYPD3-unCP-atz')
            qcdb.compare_values(-8.1985, h2kc * df['B3LYPD2-unCP-dadz']['S22-22'], digits, 'B3LYPD2-unCP-dadz')
            qcdb.compare_values(-8.1619, h2kc * df['B3LYPD3-unCP-dadz']['S22-22'], digits, 'B3LYPD3-unCP-dadz')
            qcdb.compare_values(-7.4659, h2kc * df['B3LYPD2-unCP-datz']['S22-22'], digits, 'B3LYPD2-unCP-datz')
            qcdb.compare_values(-7.4292, h2kc * df['B3LYPD3-unCP-datz']['S22-22'], digits, 'B3LYPD3-unCP-datz')
            qcdb.compare_values(-7.7993, h2kc * df['B3LYPD2-unCP-hadz']['S22-22'], digits, 'B3LYPD2-unCP-hadz')
            qcdb.compare_values(-7.7627, h2kc * df['B3LYPD3-unCP-hadz']['S22-22'], digits, 'B3LYPD3-unCP-hadz')
            qcdb.compare_values(-7.4063, h2kc * df['B3LYPD2-unCP-hatz']['S22-22'], digits, 'B3LYPD2-unCP-hatz')
            qcdb.compare_values(-7.3696, h2kc * df['B3LYPD3-unCP-hatz']['S22-22'], digits, 'B3LYPD3-unCP-hatz')
            qcdb.compare_values(-7.6638, h2kc * df['B3LYPD2-unCP-6311pg_3df_2p_']['S22-22'], digits, 'B3LYPD2-unCP-6311pg_3df_2p_')
            qcdb.compare_values(-7.6272, h2kc * df['B3LYPD3-unCP-6311pg_3df_2p_']['S22-22'], digits, 'B3LYPD3-unCP-6311pg_3df_2p_')
            qcdb.compare_values(-8.4272, h2kc * df['B3LYPD2-unCP-tz']['S22-22'], digits, 'B3LYPD2-unCP-tz')
            qcdb.compare_values(-8.3906, h2kc * df['B3LYPD3-unCP-tz']['S22-22'], digits, 'B3LYPD3-unCP-tz')
            qcdb.compare_values(-10.6481, h2kc * df['B3LYPD2-unCP-dz']['S22-22'], digits, 'B3LYPD2-unCP-dz')
            qcdb.compare_values(-10.6115, h2kc * df['B3LYPD3-unCP-dz']['S22-22'], digits, 'B3LYPD3-unCP-dz')
            qcdb.compare_values(-7.0947, h2kc * df['B3LYPXDM-unCP-adz']['S22-22'], digits, 'B3LYPXDM-unCP-adz')
            qcdb.compare_values(-6.6365, h2kc * df['B3LYPXDM-unCP-atz']['S22-22'], digits, 'B3LYPXDM-unCP-atz')

            # B2PLYP
            qcdb.compare_values(-6.0428, h2kc * df['B2PLYP-unCP-adz']['S22-22'], digits, 'B2PLYP-unCP-adz')
            qcdb.compare_values(-5.3095, h2kc * df['B2PLYP-unCP-atz']['S22-22'], digits, 'B2PLYP-unCP-atz')
            qcdb.compare_values(-8.2894, h2kc * df['B2PLYPD2-unCP-adz']['S22-22'], digits, 'B2PLYPD2-unCP-adz')
            qcdb.compare_values(-7.5561, h2kc * df['B2PLYPD2-unCP-atz']['S22-22'], digits, 'B2PLYPD2-unCP-atz')
            qcdb.compare_values(-8.2907, h2kc * df['B2PLYPD3-unCP-adz']['S22-22'], digits, 'B2PLYPD3-unCP-adz')  # new dash
            qcdb.compare_values(-7.5574, h2kc * df['B2PLYPD3-unCP-atz']['S22-22'], digits, 'B2PLYPD3-unCP-atz')  # new dash
            #qcdb.compare_values(-8.2709, h2kc * df['B2PLYPD3-unCP-adz']['S22-22'], digits, 'B2PLYPD3-unCP-adz')  # dftbench paper, old dash
            #qcdb.compare_values(-7.5376, h2kc * df['B2PLYPD3-unCP-atz']['S22-22'], digits, 'B2PLYPD3-unCP-atz')  # dftbench paper, old dash

            # B970
            qcdb.compare_values(-3.0925, h2kc * df['B970-CP-adz']['S22-22'], digits, 'B970-CP-adz')
            qcdb.compare_values(-3.1405, h2kc * df['B970-CP-atz']['S22-22'], digits, 'B970-CP-atz')
            qcdb.compare_values(-3.2000, h2kc * df['B970-CP-dz']['S22-22'], digits, 'B970-CP-dz')
            qcdb.compare_values(-3.0855, h2kc * df['B970-CP-hadz']['S22-22'], digits, 'B970-CP-hadz')
            qcdb.compare_values(-3.1304, h2kc * df['B970-CP-hatz']['S22-22'], digits, 'B970-CP-hatz')
            qcdb.compare_values(-3.1281, h2kc * df['B970-CP-tz']['S22-22'], digits, 'B970-CP-tz')
            qcdb.compare_values(-3.8891, h2kc * df['B970-unCP-adz']['S22-22'], digits, 'B970-unCP-adz')
            qcdb.compare_values(-3.3569, h2kc * df['B970-unCP-atz']['S22-22'], digits, 'B970-unCP-atz')
            qcdb.compare_values(-6.2374, h2kc * df['B970-unCP-dz']['S22-22'], digits, 'B970-unCP-dz')
            qcdb.compare_values(-3.7382, h2kc * df['B970-unCP-hadz']['S22-22'], digits, 'B970-unCP-hadz')
            qcdb.compare_values(-3.3079, h2kc * df['B970-unCP-hatz']['S22-22'], digits, 'B970-unCP-hatz')
            qcdb.compare_values(-4.2457, h2kc * df['B970-unCP-tz']['S22-22'], digits, 'B970-unCP-tz')
            qcdb.compare_values(-6.1561, h2kc * df['B970D2-CP-adz']['S22-22'], digits, 'S22-B970D2-CP-adz')
            qcdb.compare_values(-6.2041, h2kc * df['B970D2-CP-atz']['S22-22'], digits, 'S22-B970D2-CP-atz')
            qcdb.compare_values(-6.2635, h2kc * df['B970D2-CP-dz']['S22-22'], digits, 'S22-B970D2-CP-dz')
            qcdb.compare_values(-6.1491, h2kc * df['B970D2-CP-hadz']['S22-22'], digits, 'S22-B970D2-CP-hadz')
            qcdb.compare_values(-6.1939, h2kc * df['B970D2-CP-hatz']['S22-22'], digits, 'S22-B970D2-CP-hatz')
            qcdb.compare_values(-6.1917, h2kc * df['B970D2-CP-tz']['S22-22'], digits, 'S22-B970D2-CP-tz')
            qcdb.compare_values(-6.9526, h2kc * df['B970D2-unCP-adz']['S22-22'], digits, 'S22-B970D2-unCP-adz')
            qcdb.compare_values(-6.4205, h2kc * df['B970D2-unCP-atz']['S22-22'], digits, 'S22-B970D2-unCP-atz')
            qcdb.compare_values(-9.3010, h2kc * df['B970D2-unCP-dz']['S22-22'], digits, 'S22-B970D2-unCP-dz')
            qcdb.compare_values(-6.8018, h2kc * df['B970D2-unCP-hadz']['S22-22'], digits, 'S22-B970D2-unCP-hadz')
            qcdb.compare_values(-6.3715, h2kc * df['B970D2-unCP-hatz']['S22-22'], digits, 'S22-B970D2-unCP-hatz')
            qcdb.compare_values(-7.3093, h2kc * df['B970D2-unCP-tz']['S22-22'], digits, 'S22-B970D2-unCP-tz')

            # B97-D
            qcdb.compare_values(-1.6396, h2kc * df['B97-unCP-adz']['S22-22'], digits, 'B97-unCP-adz')
            qcdb.compare_values(-1.1556, h2kc * df['B97-unCP-atz']['S22-22'], digits, 'B97-unCP-atz')
            qcdb.compare_values(-1.4969, h2kc * df['B97-unCP-hadz']['S22-22'], digits, 'B97-unCP-hadz')
            qcdb.compare_values(-1.0966, h2kc * df['B97-unCP-hatz']['S22-22'], digits, 'B97-unCP-hatz')
            qcdb.compare_values(-6.7455, h2kc * df['B97D2-unCP-adz']['S22-22'], digits, 'B97D2-unCP-adz')
            qcdb.compare_values(-7.0599, h2kc * df['B97D3-unCP-adz']['S22-22'], digits, 'B97D3-unCP-adz')
            qcdb.compare_values(-6.2616, h2kc * df['B97D2-unCP-atz']['S22-22'], digits, 'B97D2-unCP-atz')
            qcdb.compare_values(-6.5760, h2kc * df['B97D3-unCP-atz']['S22-22'], digits, 'B97D3-unCP-atz')
            qcdb.compare_values(-6.9173, h2kc * df['B97D3-unCP-hadz']['S22-22'], digits, 'B97D3-unCP-hadz')
            qcdb.compare_values(-6.5170, h2kc * df['B97D3-unCP-hatz']['S22-22'], digits, 'B97D3-unCP-hatz')

            # BP86
            qcdb.compare_values(-2.2746, h2kc * df['BP86-CP-adz']['S22-22'], digits, 'BP86-CP-adz')
            qcdb.compare_values(-2.3338, h2kc * df['BP86-CP-atz']['S22-22'], digits, 'BP86-CP-atz')
            qcdb.compare_values(-2.3254, h2kc * df['BP86-CP-dz']['S22-22'], digits, 'BP86-CP-dz')
            qcdb.compare_values(-2.2529, h2kc * df['BP86-CP-hadz']['S22-22'], digits, 'BP86-CP-hadz')
            qcdb.compare_values(-2.3292, h2kc * df['BP86-CP-hatz']['S22-22'], digits, 'BP86-CP-hatz')
            qcdb.compare_values(-2.3577, h2kc * df['BP86-CP-tz']['S22-22'], digits, 'BP86-CP-tz')
            qcdb.compare_values(-3.0384, h2kc * df['BP86-unCP-adz']['S22-22'], digits, 'BP86-unCP-adz')
            qcdb.compare_values(-2.6030, h2kc * df['BP86-unCP-atz']['S22-22'], digits, 'BP86-unCP-atz')
            qcdb.compare_values(-5.4818, h2kc * df['BP86-unCP-dz']['S22-22'], digits, 'BP86-unCP-dz')
            qcdb.compare_values(-2.8822, h2kc * df['BP86-unCP-hadz']['S22-22'], digits, 'BP86-unCP-hadz')
            qcdb.compare_values(-2.5505, h2kc * df['BP86-unCP-hatz']['S22-22'], digits, 'BP86-unCP-hatz')
            qcdb.compare_values(-3.4885, h2kc * df['BP86-unCP-tz']['S22-22'], digits, 'BP86-unCP-tz')
            qcdb.compare_values(-6.5636, h2kc * df['BP86D2-CP-adz']['S22-22'], digits, 'BP86D2-CP-adz')
            qcdb.compare_values(-6.6228, h2kc * df['BP86D2-CP-atz']['S22-22'], digits, 'BP86D2-CP-atz')
            qcdb.compare_values(-6.6145, h2kc * df['BP86D2-CP-dz']['S22-22'], digits, 'BP86D2-CP-dz')
            qcdb.compare_values(-6.5419, h2kc * df['BP86D2-CP-hadz']['S22-22'], digits, 'BP86D2-CP-hadz')
            qcdb.compare_values(-6.6182, h2kc * df['BP86D2-CP-hatz']['S22-22'], digits, 'BP86D2-CP-hatz')
            qcdb.compare_values(-6.6467, h2kc * df['BP86D2-CP-tz']['S22-22'], digits, 'BP86D2-CP-tz')
            qcdb.compare_values(-7.3274, h2kc * df['BP86D2-unCP-adz']['S22-22'], digits, 'BP86D2-unCP-adz')
            qcdb.compare_values(-7.9997, h2kc * df['BP86D3-unCP-adz']['S22-22'], digits, 'BP86D3-unCP-adz')
            qcdb.compare_values(-6.8920, h2kc * df['BP86D2-unCP-atz']['S22-22'], digits, 'BP86D2-unCP-atz')
            qcdb.compare_values(-7.5644, h2kc * df['BP86D3-unCP-atz']['S22-22'], digits, 'BP86D3-unCP-atz')
            qcdb.compare_values(-9.7708, h2kc * df['BP86D2-unCP-dz']['S22-22'], digits, 'BP86D2-unCP-dz')
            qcdb.compare_values(-7.1712, h2kc * df['BP86D2-unCP-hadz']['S22-22'], digits, 'BP86D2-unCP-hadz')
            qcdb.compare_values(-6.8395, h2kc * df['BP86D2-unCP-hatz']['S22-22'], digits, 'BP86D2-unCP-hatz')
            qcdb.compare_values(-7.7775, h2kc * df['BP86D2-unCP-tz']['S22-22'], digits, 'BP86D2-unCP-tz')

            # wB97X-D
            qcdb.compare_values(-6.9661, h2kc * df['WB97XD-CP-adz']['S22-22'], digits, 'WB97XD-CP-adz')
            qcdb.compare_values(-6.9069, h2kc * df['WB97XD-CP-atz']['S22-22'], digits, 'WB97XD-CP-atz')
            qcdb.compare_values(-6.8996, h2kc * df['WB97XD-CP-dz']['S22-22'], digits, 'WB97XD-CP-dz')
            qcdb.compare_values(-6.9392, h2kc * df['WB97XD-CP-hadz']['S22-22'], digits, 'WB97XD-CP-hadz')
            qcdb.compare_values(-6.9245, h2kc * df['WB97XD-CP-hatz']['S22-22'], digits, 'WB97XD-CP-hatz')
            qcdb.compare_values(-6.9123, h2kc * df['WB97XD-CP-tz']['S22-22'], digits, 'WB97XD-CP-tz')
            qcdb.compare_values(-7.7283, h2kc * df['WB97XD-unCP-adz']['S22-22'], digits, 'WB97XD-unCP-adz')
            qcdb.compare_values(-7.2139, h2kc * df['WB97XD-unCP-atz']['S22-22'], digits, 'WB97XD-unCP-atz')
            qcdb.compare_values(-9.6968, h2kc * df['WB97XD-unCP-dz']['S22-22'], digits, 'WB97XD-unCP-dz')
            qcdb.compare_values(-7.5543, h2kc * df['WB97XD-unCP-hadz']['S22-22'], digits, 'WB97XD-unCP-hadz')
            qcdb.compare_values(-7.1715, h2kc * df['WB97XD-unCP-hatz']['S22-22'], digits, 'WB97XD-unCP-hatz')
            qcdb.compare_values(-7.9090, h2kc * df['WB97XD-unCP-tz']['S22-22'], digits, 'WB97XD-unCP-tz')

            # M05-2X
            qcdb.compare_values(-5.9556, h2kc * df['M052X-CP-adz']['S22-22'], digits, 'M052X-CP-adz')
            qcdb.compare_values(-6.0724, h2kc * df['M052X-CP-atz']['S22-22'], digits, 'M052X-CP-atz')
            qcdb.compare_values(-5.8742, h2kc * df['M052X-CP-dz']['S22-22'], digits, 'M052X-CP-dz')
            qcdb.compare_values(-5.9155, h2kc * df['M052X-CP-hadz']['S22-22'], digits, 'M052X-CP-hadz')
            qcdb.compare_values(-6.0799, h2kc * df['M052X-CP-hatz']['S22-22'], digits, 'M052X-CP-hatz')
            qcdb.compare_values(-6.0215, h2kc * df['M052X-CP-tz']['S22-22'], digits, 'M052X-CP-tz')
            qcdb.compare_values(-6.7954, h2kc * df['M052X-unCP-adz']['S22-22'], digits, 'M052X-unCP-adz')
            qcdb.compare_values(-6.4217, h2kc * df['M052X-unCP-atz']['S22-22'], digits, 'M052X-unCP-atz')
            qcdb.compare_values(-8.6857, h2kc * df['M052X-unCP-dz']['S22-22'], digits, 'M052X-unCP-dz')
            qcdb.compare_values(-6.5930, h2kc * df['M052X-unCP-hadz']['S22-22'], digits, 'M052X-unCP-hadz')
            qcdb.compare_values(-6.3553, h2kc * df['M052X-unCP-hatz']['S22-22'], digits, 'M052X-unCP-hatz')
            qcdb.compare_values(-6.9496, h2kc * df['M052X-unCP-tz']['S22-22'], digits, 'M052X-unCP-tz')
            qcdb.compare_values(-7.6826, h2kc * df['M052XD3-unCP-adz']['S22-22'], digits, 'M052XD3-unCP-adz')
            qcdb.compare_values(-7.3088, h2kc * df['M052XD3-unCP-atz']['S22-22'], digits, 'M052XD3-unCP-atz')

            # M06-2X
            qcdb.compare_values(-6.5146, h2kc * df['M062X-CP-adz']['S22-22'], digits, 'M062X-CP-adz')
            qcdb.compare_values(-6.5581, h2kc * df['M062X-CP-atz']['S22-22'], digits, 'M062X-CP-atz')
            qcdb.compare_values(-6.3411, h2kc * df['M062X-CP-dz']['S22-22'], digits, 'M062X-CP-dz')
            qcdb.compare_values(-6.4566, h2kc * df['M062X-CP-hadz']['S22-22'], digits, 'M062X-CP-hadz')
            qcdb.compare_values(-6.5696, h2kc * df['M062X-CP-hatz']['S22-22'], digits, 'M062X-CP-hatz')
            qcdb.compare_values(-6.4499, h2kc * df['M062X-CP-tz']['S22-22'], digits, 'M062X-CP-tz')
            qcdb.compare_values(-7.3393, h2kc * df['M062X-unCP-adz']['S22-22'], digits, 'M062X-unCP-adz')
            qcdb.compare_values(-6.8876, h2kc * df['M062X-unCP-atz']['S22-22'], digits, 'M062X-unCP-atz')
            qcdb.compare_values(-9.1232, h2kc * df['M062X-unCP-dz']['S22-22'], digits, 'M062X-unCP-dz')
            qcdb.compare_values(-7.1083, h2kc * df['M062X-unCP-hadz']['S22-22'], digits, 'M062X-unCP-hadz')
            qcdb.compare_values(-6.8307, h2kc * df['M062X-unCP-hatz']['S22-22'], digits, 'M062X-unCP-hatz')
            qcdb.compare_values(-7.4349, h2kc * df['M062X-unCP-tz']['S22-22'], digits, 'M062X-unCP-tz')
            qcdb.compare_values(-7.8457, h2kc * df['M062XD3-unCP-adz']['S22-22'], digits, 'M062XD3-unCP-adz')
            qcdb.compare_values(-7.3941, h2kc * df['M062XD3-unCP-atz']['S22-22'], digits, 'M062XD3-unCP-atz')

            # PBE
            qcdb.compare_values(-3.8225, h2kc * df['PBE-CP-adz']['S22-22'], digits, 'PBE-CP-adz')
            qcdb.compare_values(-3.8524, h2kc * df['PBE-CP-atz']['S22-22'], digits, 'PBE-CP-atz')
            qcdb.compare_values(-3.8816, h2kc * df['PBE-CP-dz']['S22-22'], digits, 'PBE-CP-dz')
            qcdb.compare_values(-3.8091, h2kc * df['PBE-CP-hadz']['S22-22'], digits, 'PBE-CP-hadz')
            qcdb.compare_values(-3.8445, h2kc * df['PBE-CP-hatz']['S22-22'], digits, 'PBE-CP-hatz')
            qcdb.compare_values(-3.8422, h2kc * df['PBE-CP-tz']['S22-22'], digits, 'PBE-CP-tz')
            qcdb.compare_values(-4.5854, h2kc * df['PBE-unCP-adz']['S22-22'], digits, 'PBE-unCP-adz')
            qcdb.compare_values(-4.0802, h2kc * df['PBE-unCP-atz']['S22-22'], digits, 'PBE-unCP-atz')
            qcdb.compare_values(-7.3095, h2kc * df['PBE-unCP-dz']['S22-22'], digits, 'PBE-unCP-dz')
            qcdb.compare_values(-4.4503, h2kc * df['PBE-unCP-hadz']['S22-22'], digits, 'PBE-unCP-hadz')
            qcdb.compare_values(-4.0317, h2kc * df['PBE-unCP-hatz']['S22-22'], digits, 'PBE-unCP-hatz')
            qcdb.compare_values(-5.1308, h2kc * df['PBE-unCP-tz']['S22-22'], digits, 'PBE-unCP-tz')
            qcdb.compare_values(-6.8861, h2kc * df['PBED2-CP-adz']['S22-22'], digits, 'PBED2-CP-adz')
            qcdb.compare_values(-6.6763, h2kc * df['PBED3-CP-adz']['S22-22'], digits, 'PBED3-CP-adz')
            qcdb.compare_values(-6.9160, h2kc * df['PBED2-CP-atz']['S22-22'], digits, 'PBED2-CP-atz')
            qcdb.compare_values(-6.7063, h2kc * df['PBED3-CP-atz']['S22-22'], digits, 'PBED3-CP-atz')
            qcdb.compare_values(-6.9451, h2kc * df['PBED2-CP-dz']['S22-22'], digits, 'PBED2-CP-dz')
            qcdb.compare_values(-6.7354, h2kc * df['PBED3-CP-dz']['S22-22'], digits, 'PBED3-CP-dz')
            qcdb.compare_values(-6.8727, h2kc * df['PBED2-CP-hadz']['S22-22'], digits, 'PBED2-CP-hadz')
            qcdb.compare_values(-6.6630, h2kc * df['PBED3-CP-hadz']['S22-22'], digits, 'PBED3-CP-hadz')
            qcdb.compare_values(-6.9081, h2kc * df['PBED2-CP-hatz']['S22-22'], digits, 'PBED2-CP-hatz')
            qcdb.compare_values(-6.6983, h2kc * df['PBED3-CP-hatz']['S22-22'], digits, 'PBED3-CP-hatz')
            qcdb.compare_values(-6.9057, h2kc * df['PBED2-CP-tz']['S22-22'], digits, 'PBED2-CP-tz')
            qcdb.compare_values(-6.6960, h2kc * df['PBED3-CP-tz']['S22-22'], digits, 'PBED3-CP-tz')
            qcdb.compare_values(-7.6490, h2kc * df['PBED2-unCP-adz']['S22-22'], digits, 'PBED2-unCP-adz')
            qcdb.compare_values(-7.4392, h2kc * df['PBED3-unCP-adz']['S22-22'], digits, 'PBED3-unCP-adz')
            qcdb.compare_values(-7.1438, h2kc * df['PBED2-unCP-atz']['S22-22'], digits, 'PBED2-unCP-atz')
            qcdb.compare_values(-6.9340, h2kc * df['PBED3-unCP-atz']['S22-22'], digits, 'PBED3-unCP-atz')
            qcdb.compare_values(-7.5139, h2kc * df['PBED2-unCP-hadz']['S22-22'], digits, 'PBED2-unCP-hadz')
            qcdb.compare_values(-7.3041, h2kc * df['PBED3-unCP-hadz']['S22-22'], digits, 'PBED3-unCP-hadz')
            qcdb.compare_values(-7.0952, h2kc * df['PBED2-unCP-hatz']['S22-22'], digits, 'PBED2-unCP-hatz')
            qcdb.compare_values(-6.8855, h2kc * df['PBED3-unCP-hatz']['S22-22'], digits, 'PBED3-unCP-hatz')
            qcdb.compare_values(-8.1944, h2kc * df['PBED2-unCP-tz']['S22-22'], digits, 'PBED2-unCP-tz')
            qcdb.compare_values(-7.9847, h2kc * df['PBED3-unCP-tz']['S22-22'], digits, 'PBED3-unCP-tz')
            qcdb.compare_values(-10.3731, h2kc * df['PBED2-unCP-dz']['S22-22'], digits, 'PBED2-unCP-dz')
            qcdb.compare_values(-10.1634, h2kc * df['PBED3-unCP-dz']['S22-22'], digits, 'PBED3-unCP-dz')

            # PBE0
            qcdb.compare_values(-4.9388, h2kc * df['PBE0-unCP-adz']['S22-22'], digits, 'PBE0-unCP-adz')
            qcdb.compare_values(-4.3888, h2kc * df['PBE0-unCP-atz']['S22-22'], digits, 'PBE0-unCP-atz')
            qcdb.compare_values(-7.3897, h2kc * df['PBE0D2-unCP-adz']['S22-22'], digits, 'PBE0D2-unCP-adz')
            qcdb.compare_values(-6.8397, h2kc * df['PBE0D2-unCP-atz']['S22-22'], digits, 'PBE0D2-unCP-atz')
            qcdb.compare_values(-7.7887, h2kc * df['PBE0D3-unCP-adz']['S22-22'], digits, 'PBE0D3-unCP-adz')  # new dash
            qcdb.compare_values(-7.2387, h2kc * df['PBE0D3-unCP-atz']['S22-22'], digits, 'PBE0D3-unCP-atz')  # new dash
            #qcdb.compare_values(-7.8275, h2kc * df['PBE0D3-unCP-adz']['S22-22'], digits, 'PBE0D3-unCP-adz')  # dftbench paper, old dash
            #qcdb.compare_values(-7.2775, h2kc * df['PBE0D3-unCP-atz']['S22-22'], digits, 'PBE0D3-unCP-atz')  # dftbench paper, old dash

            # XYG3
            qcdb.compare_values(-7.2745, h2kc * df['XYG3-unCP-adz']['S22-22'], digits, 'XYG3-unCP-adz')
            qcdb.compare_values(-6.4390, h2kc * df['XYG3-unCP-atz']['S22-22'], digits, 'XYG3-unCP-atz')
            qcdb.compare_values(-7.5166, h2kc * df['XYG3-unCP-6311pg_3df_2p_']['S22-22'], digits, 'XYG3-unCP-6311pg_3df_2p_')
            qcdb.compare_values(-6.7097, h2kc * df['XYG3-unCP-6311ppg_3df_2p_']['S22-22'], digits, 'XYG3-unCP-6311ppg_3df_2p_')

        elif project == 'parenq':
            print """\n        Checks from parenq project:\n"""
            digits = 4

            qcdb.compare_values(-3.7422, h2kc * df['HF-CP-atz']['A24-4'], digits, 'HF-CP-atz')
            qcdb.compare_values(-3.7558, h2kc * df['HF-CP-adz']['A24-4'], digits, 'HF-CP-adz')
            qcdb.compare_values(-3.8121, h2kc * df['HF-CP-hadz']['A24-4'], digits, 'HF-CP-hadz')
            qcdb.compare_values(-3.7479, h2kc * df['HF-CP-jadz']['A24-4'], digits, 'HF-CP-jadz')

            qcdb.compare_values(-4.2839, h2kc * df['CCSDT-fno1e5-CP-atz']['A24-4'], digits, 'CCSDT-fno1e5-CP-atz')
            qcdb.compare_values(-3.8546, h2kc * df['CCSDT-fno1e4-CP-atz']['A24-4'], digits, 'CCSDT-fno1e4-CP-atz')
            qcdb.compare_values(-3.9771, h2kc * df['CCSDT-full-CP-adz']['A24-4'], digits, 'CCSDT-full-CP-adz')
            qcdb.compare_values(-3.9894, h2kc * df['CCSDT-fno1e6-CP-adz']['A24-4'], digits, 'CCSDT-fno1e6-CP-adz')
            qcdb.compare_values(-4.0386, h2kc * df['CCSDT-fno1e5-CP-adz']['A24-4'], digits, 'CCSDT-fno1e5-CP-adz')
            qcdb.compare_values(-4.1458, h2kc * df['CCSDT-fno5e5-CP-adz']['A24-4'], digits, 'CCSDT-fno5e5-CP-adz')
            qcdb.compare_values(-3.9510, h2kc * df['CCSDT-fno1e4-CP-adz']['A24-4'], digits, 'CCSDT-fno1e4-CP-adz')
            qcdb.compare_values(-3.9510, h2kc * df['CCSDT-fno1e4-mrcc-CP-adz']['A24-4'], digits, 'CCSDT-fno1e4-mrcc-CP-adz')
            qcdb.compare_values(-3.8526, h2kc * df['CCSDT-full-CP-hadz']['A24-4'], digits, 'CCSDT-full-CP-hadz')
            qcdb.compare_values(-3.4657, h2kc * df['CCSDT-full-CP-jadz']['A24-4'], digits, 'CCSDT-full-CP-jadz')

            qcdb.compare_values(-4.2905, h2kc * df['CCSDTQ-fno1e5-CP-atz']['A24-4'], digits, 'CCSDTQ-fno1e5-CP-atz')
            qcdb.compare_values(-3.8583, h2kc * df['CCSDTQ-fno1e4-CP-atz']['A24-4'], digits, 'CCSDTQ-fno1e4-CP-atz')
            qcdb.compare_values(-3.9703, h2kc * df['CCSDTQ-full-CP-adz']['A24-4'], digits, 'CCSDTQ-full-CP-adz')
            qcdb.compare_values(-3.9836, h2kc * df['CCSDTQ-fno1e6-CP-adz']['A24-4'], digits, 'CCSDTQ-fno1e6-CP-adz')
            qcdb.compare_values(-4.0340, h2kc * df['CCSDTQ-fno1e5-CP-adz']['A24-4'], digits, 'CCSDTQ-fno1e5-CP-adz')
            qcdb.compare_values(-4.1418, h2kc * df['CCSDTQ-fno5e5-CP-adz']['A24-4'], digits, 'CCSDTQ-fno5e5-CP-adz')
            qcdb.compare_values(-3.9456, h2kc * df['CCSDTQ-fno1e4-CP-adz']['A24-4'], digits, 'CCSDTQ-fno1e4-CP-adz')
            qcdb.compare_values(-3.9964, h2kc * df['CCSDTQ-full-CP-hadz']['A24-4'], digits, 'CCSDTQ-full-CP-hadz')
            qcdb.compare_values(-3.4632, h2kc * df['CCSDTQ-full-CP-jadz']['A24-4'], digits, 'CCSDTQ-full-CP-jadz')

            qcdb.compare_values(-3.0433, h2kc * df['CCSDTQ-fno1e4-CP-adtz']['A24-9'], digits, 'CCSDTQ-fno1e4-CP-adtz')

        elif project == 'f12dilabio':
            print """\n        Checks from f12dilabio project:\n"""
            digits = 4

            qcdb.compare_values(-2.4749, h2kc * df['HFCABS-CP-adz']['A24-9'], digits, 'HFCABS-CP-adz')
            qcdb.compare_values(-4.5972, h2kc * df['CCSDTAF12-CP-adz']['A24-9'], digits, 'CCSDTAF12-CP-adz')
            qcdb.compare_values(-4.4354, h2kc * df['CCSDTBF12-CP-adz']['A24-9'], digits, 'CCSDTBF12-CP-adz')
            qcdb.compare_values(-4.4726, h2kc * df['CCSDTCF12-CP-adz']['A24-9'], digits, 'CCSDTCF12-CP-adz')  # corr, translate error b/c tz scf
            qcdb.compare_values(-4.5830, h2kc * df['DWCCSDTF12-CP-adz']['A24-9'], digits, 'DWCCSDTF12-CP-adz')
            qcdb.compare_values(-4.5913, h2kc * df['DWCCSDTF12-CP-atz']['A24-9'], digits, 'DWCCSDTF12-CP-atz')  # corr, not separate atz def in reapsets
            qcdb.compare_values(-4.5669, h2kc * df['MP2F12-CP-adz']['A24-9'], digits, 'MP2F12-CP-adz')
            qcdb.compare_values(-4.5767, h2kc * df['CCSDTAF12-CP-a5z']['A24-9'], digits, 'CCSDTAF12-CP-a5z')
            qcdb.compare_values(-4.5855, h2kc * df['CCSDTAF12-CP-aqz']['A24-9'], digits, 'CCSDTAF12-CP-aqz')
            qcdb.compare_values(-4.6042, h2kc * df['CCSDTAF12-CP-atz']['A24-9'], digits, 'CCSDTAF12-CP-atz')
            qcdb.compare_values(-4.3559, h2kc * df['CCSDTAF12-CP-dzf12']['A24-9'], digits, 'CCSDTAF12-CP-dzf12')  # corr, translate error b/c tz scf
            qcdb.compare_values(-4.2422, h2kc * df['CCSDTBF12-CP-dzf12']['A24-9'], digits, 'CCSDTBF12-CP-dzf12')  # corr, translate error b/c tz scf
            qcdb.compare_values(-4.2744, h2kc * df['CCSDTCF12-CP-dzf12']['A24-9'], digits, 'CCSDTCF12-CP-dzf12')  # corr, translate error b/c tz scf
            qcdb.compare_values(-2.4582, h2kc * df['HFCABS-CP-dzf12']['A24-9'], digits, 'HFCABS-CP-dzf12')  # added
            #qcdb.compare_values(, h2kc * df['HFCABS-CP-tzf12']['A24-9'], digits, 'HFCABS-CP-tzf12')  # added
            #qcdb.compare_values(-2.4692, h2kc * df['HFCABS-CP-qzf12']['A24-9'], digits, 'HFCABS-CP-qzf12')  # added
            qcdb.compare_values(-4.5204, h2kc * df['CCSDTAF12-CP-tzf12']['A24-9'], digits, 'CCSDTAF12-CP-tzf12')
            qcdb.compare_values(-4.5692, h2kc * df['CCSDTAF12-CP-qzf12']['A24-9'], digits, 'CCSDTAF12-CP-qzf12')
            qcdb.compare_values(-4.5556, h2kc * df['CCSDTBF12-CP-a5z']['A24-9'], digits, 'CCSDTBF12-CP-a5z')
            qcdb.compare_values(-4.5495, h2kc * df['CCSDTBF12-CP-aqz']['A24-9'], digits, 'CCSDTBF12-CP-aqz')
            qcdb.compare_values(-4.5365, h2kc * df['CCSDTBF12-CP-atz']['A24-9'], digits, 'CCSDTBF12-CP-atz')
            qcdb.compare_values(-4.6084, h2kc * df['CCSDTAF12-CP-adtz']['A24-9'], digits, 'CCSDTAF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5673, h2kc * df['CCSDTAF12-CP-aq5z']['A24-9'], digits, 'CCSDTAF12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5705, h2kc * df['CCSDTAF12-CP-atqz']['A24-9'], digits, 'CCSDTAF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5851, h2kc * df['CCSDTAF12-CP-dtzf12']['A24-9'], digits, 'CCSDTAF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.6017, h2kc * df['CCSDTAF12-CP-tqzf12']['A24-9'], digits, 'CCSDTAF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5804, h2kc * df['CCSDTBF12-CP-adtz']['A24-9'], digits, 'CCSDTBF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5619, h2kc * df['CCSDTBF12-CP-aq5z']['A24-9'], digits, 'CCSDTBF12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5576, h2kc * df['CCSDTBF12-CP-atqz']['A24-9'], digits, 'CCSDTBF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5547, h2kc * df['CCSDTBF12-CP-dtzf12']['A24-9'], digits, 'CCSDTBF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5873, h2kc * df['CCSDTBF12-CP-tqzf12']['A24-9'], digits, 'CCSDTBF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5752, h2kc * df['CCSDTCF12-CP-adtz']['A24-9'], digits, 'CCSDTCF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5447, h2kc * df['CCSDTCF12-CP-atqz']['A24-9'], digits, 'CCSDTCF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5462, h2kc * df['CCSDTCF12-CP-dtzf12']['A24-9'], digits, 'CCSDTCF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5756, h2kc * df['CCSDTCF12-CP-tqzf12']['A24-9'], digits, 'CCSDTCF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5377, h2kc * df['CCSDTBF12-CP-qzf12']['A24-9'], digits, 'CCSDTBF12-CP-qzf12')
            qcdb.compare_values(-4.4653, h2kc * df['CCSDTBF12-CP-tzf12']['A24-9'], digits, 'CCSDTBF12-CP-tzf12')
            qcdb.compare_values(-4.5451, h2kc * df['CCSDTCF12-CP-aqz']['A24-9'], digits, 'CCSDTCF12-CP-aqz')
            qcdb.compare_values(-4.5439, h2kc * df['CCSDTCF12-CP-atz']['A24-9'], digits, 'CCSDTCF12-CP-atz')
            qcdb.compare_values(-4.5324, h2kc * df['CCSDTCF12-CP-qzf12']['A24-9'], digits, 'CCSDTCF12-CP-qzf12')
            qcdb.compare_values(-4.4688, h2kc * df['CCSDTCF12-CP-tzf12']['A24-9'], digits, 'CCSDTCF12-CP-tzf12')
            qcdb.compare_values(-2.4739, h2kc * df['HFCABS-CP-a5z']['A24-9'], digits, 'HFCABS-CP-a5z')
            qcdb.compare_values(-2.4737, h2kc * df['HFCABS-CP-aqz']['A24-9'], digits, 'HFCABS-CP-aqz')
            qcdb.compare_values(-2.4719, h2kc * df['HFCABS-CP-atz']['A24-9'], digits, 'HFCABS-CP-atz')
            qcdb.compare_values(-4.6070, h2kc * df['MP2F12-CP-atz']['A24-9'], digits, 'MP2F12-CP-atz')
            qcdb.compare_values(-4.5927, h2kc * df['MP2F12-CP-aqz']['A24-9'], digits, 'MP2F12-CP-aqz')
            qcdb.compare_values(-4.5954, h2kc * df['MP2F12-CP-a5z']['A24-9'], digits, 'MP2F12-CP-a5z')
            qcdb.compare_values(-4.6251, h2kc * df['MP2F12-CP-adtz']['A24-9'], digits, 'MP2F12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5980, h2kc * df['MP2F12-CP-aq5z']['A24-9'], digits, 'MP2F12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.5810, h2kc * df['MP2F12-CP-atqz']['A24-9'], digits, 'MP2F12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
            qcdb.compare_values(-4.6099, h2kc * df['CCSDTAF12-CP-hill1_adtz']['A24-9'], digits, 'CCSDTAF12-CP-hill1_adtz')
            qcdb.compare_values(-4.5769, h2kc * df['CCSDTAF12-CP-hill1_atqz']['A24-9'], digits, 'CCSDTAF12-CP-hill1_atqz')
            qcdb.compare_values(-4.5722, h2kc * df['CCSDTAF12-CP-hill1_aq5z']['A24-9'], digits, 'CCSDTAF12-CP-hill1_aq5z')
            qcdb.compare_values(-4.5800, h2kc * df['CCSDTAF12-CP-hill1_dtzf12']['A24-9'], digits, 'CCSDTAF12-CP-hill1_dtzf12')
            qcdb.compare_values(-4.5854, h2kc * df['CCSDTAF12-CP-hill1_tqzf12']['A24-9'], digits, 'CCSDTAF12-CP-hill1_tqzf12')
            qcdb.compare_values(-4.5965, h2kc * df['CCSDTBF12-CP-hill1_adtz']['A24-9'], digits, 'CCSDTBF12-CP-hill1_adtz')
            qcdb.compare_values(-4.5541, h2kc * df['CCSDTBF12-CP-hill1_atqz']['A24-9'], digits, 'CCSDTBF12-CP-hill1_atqz')
            qcdb.compare_values(-4.5586, h2kc * df['CCSDTBF12-CP-hill1_aq5z']['A24-9'], digits, 'CCSDTBF12-CP-hill1_aq5z')
            qcdb.compare_values(-4.5476, h2kc * df['CCSDTBF12-CP-hill1_dtzf12']['A24-9'], digits, 'CCSDTBF12-CP-hill1_dtzf12')
            qcdb.compare_values(-4.5624, h2kc * df['CCSDTBF12-CP-hill1_tqzf12']['A24-9'], digits, 'CCSDTBF12-CP-hill1_tqzf12')
            qcdb.compare_values(-4.5867, h2kc * df['CCSDTCF12-CP-hill1_adtz']['A24-9'], digits, 'CCSDTCF12-CP-hill1_adtz')
            qcdb.compare_values(-4.5449, h2kc * df['CCSDTCF12-CP-hill1_atqz']['A24-9'], digits, 'CCSDTCF12-CP-hill1_atqz')
            qcdb.compare_values(-4.5400, h2kc * df['CCSDTCF12-CP-hill1_dtzf12']['A24-9'], digits, 'CCSDTCF12-CP-hill1_dtzf12')
            qcdb.compare_values(-4.5539, h2kc * df['CCSDTCF12-CP-hill1_tqzf12']['A24-9'], digits, 'CCSDTCF12-CP-hill1_tqzf12')
            qcdb.compare_values(-4.1343, h2kc * df['CCSDAF12-CP-atz']['A24-9'], digits, 'CCSDAF12-CP-atz')
            qcdb.compare_values(-4.1111, h2kc * df['CCSDAF12-CP-a5z']['A24-9'], digits, 'CCSDAF12-CP-a5z')
            qcdb.compare_values(-4.1183, h2kc * df['CCSDAF12-CP-aqz']['A24-9'], digits, 'CCSDAF12-CP-aqz')
            qcdb.compare_values(-4.5610, h2kc * df['CCSDTNSAF12-CP-atz']['A24-9'], digits, 'CCSDTNSAF12-CP-atz')
            qcdb.compare_values(-4.4933, h2kc * df['CCSDTNSBF12-CP-atz']['A24-9'], digits, 'CCSDTNSBF12-CP-atz')
            qcdb.compare_values(-4.5892, h2kc * df['CCSDTNSBF12-CP-hill2_adtz']['A24-9'], digits, 'CCSDTNSBF12-CP-hill2_adtz')

        elif project == 'dilabio':
            print """\n        Checks from dilabio project:\n"""
            digits = 4

            qcdb.compare_values(-4.4367, h2kc * df['CCSD-CP-a56z']['A24-4'], digits, 'CCSD-CP-a56z')
            qcdb.compare_values(-4.3994, h2kc * df['CCSD-CP-a5z']['A24-4'], digits, 'CCSD-CP-a5z')
            qcdb.compare_values(-4.4157, h2kc * df['CCSD-CP-a6z']['A24-4'], digits, 'CCSD-CP-a6z')
            qcdb.compare_values(-4.2873, h2kc * df['CCSD-CP-adtz']['A24-4'], digits, 'CCSD-CP-adtz')
            qcdb.compare_values(-3.8727, h2kc * df['CCSD-CP-adz']['A24-4'], digits, 'CCSD-CP-adz')
            qcdb.compare_values(-4.4372, h2kc * df['CCSD-CP-aq5z']['A24-4'], digits, 'CCSD-CP-aq5z')
            qcdb.compare_values(-4.3586, h2kc * df['CCSD-CP-aqz']['A24-4'], digits, 'CCSD-CP-aqz')
            qcdb.compare_values(-4.4495, h2kc * df['CCSD-CP-atqz']['A24-4'], digits, 'CCSD-CP-atqz')
            qcdb.compare_values(-4.1604, h2kc * df['CCSD-CP-atz']['A24-4'], digits, 'CCSD-CP-atz')
            qcdb.compare_values(-4.4246, h2kc * df['CCSD-unCP-a56z']['A24-4'], digits, 'CCSD-unCP-a56z')
            qcdb.compare_values(-4.5158, h2kc * df['CCSD-unCP-a5z']['A24-4'], digits, 'CCSD-unCP-a5z')
            qcdb.compare_values(-4.4739, h2kc * df['CCSD-unCP-a6z']['A24-4'], digits, 'CCSD-unCP-a6z')
            qcdb.compare_values(-4.7420, h2kc * df['CCSD-unCP-adtz']['A24-4'], digits, 'CCSD-unCP-adtz')
            qcdb.compare_values(-4.5568, h2kc * df['CCSD-unCP-adz']['A24-4'], digits, 'CCSD-unCP-adz')
            qcdb.compare_values(-4.4897, h2kc * df['CCSD-unCP-aq5z']['A24-4'], digits, 'CCSD-unCP-aq5z')
            qcdb.compare_values(-4.5825, h2kc * df['CCSD-unCP-aqz']['A24-4'], digits, 'CCSD-unCP-aqz')
            qcdb.compare_values(-4.5278, h2kc * df['CCSD-unCP-atqz']['A24-4'], digits, 'CCSD-unCP-atqz')
            qcdb.compare_values(-4.6469, h2kc * df['CCSD-unCP-atz']['A24-4'], digits, 'CCSD-unCP-atz')
            qcdb.compare_values(-4.5920, h2kc * df['CCSDT-CP-a56z']['A24-4'], digits, 'CCSDT-CP-a56z')
            qcdb.compare_values(-4.5521, h2kc * df['CCSDT-CP-a5z']['A24-4'], digits, 'CCSDT-CP-a5z')
            qcdb.compare_values(-4.5695, h2kc * df['CCSDT-CP-a6z']['A24-4'], digits, 'CCSDT-CP-a6z')
            qcdb.compare_values(-4.4422, h2kc * df['CCSDT-CP-adtz']['A24-4'], digits, 'CCSDT-CP-adtz')
            qcdb.compare_values(-3.9771, h2kc * df['CCSDT-CP-adz']['A24-4'], digits, 'CCSDT-CP-adz')
            qcdb.compare_values(-4.5928, h2kc * df['CCSDT-CP-aq5z']['A24-4'], digits, 'CCSDT-CP-aq5z')
            qcdb.compare_values(-4.5086, h2kc * df['CCSDT-CP-aqz']['A24-4'], digits, 'CCSDT-CP-aqz')
            qcdb.compare_values(-4.6067, h2kc * df['CCSDT-CP-atqz']['A24-4'], digits, 'CCSDT-CP-atqz')
            qcdb.compare_values(-4.3004, h2kc * df['CCSDT-CP-atz']['A24-4'], digits, 'CCSDT-CP-atz')
            qcdb.compare_values(-4.5785, h2kc * df['CCSDT-unCP-a56z']['A24-4'], digits, 'CCSDT-unCP-a56z')
            qcdb.compare_values(-4.6756, h2kc * df['CCSDT-unCP-a5z']['A24-4'], digits, 'CCSDT-unCP-a5z')
            qcdb.compare_values(-4.6312, h2kc * df['CCSDT-unCP-a6z']['A24-4'], digits, 'CCSDT-unCP-a6z')
            qcdb.compare_values(-4.8890, h2kc * df['CCSDT-unCP-adtz']['A24-4'], digits, 'CCSDT-unCP-adtz')
            qcdb.compare_values(-4.7587, h2kc * df['CCSDT-unCP-adz']['A24-4'], digits, 'CCSDT-unCP-adz')
            qcdb.compare_values(-4.6504, h2kc * df['CCSDT-unCP-aq5z']['A24-4'], digits, 'CCSDT-unCP-aq5z')
            qcdb.compare_values(-4.7415, h2kc * df['CCSDT-unCP-aqz']['A24-4'], digits, 'CCSDT-unCP-aqz')
            qcdb.compare_values(-4.6836, h2kc * df['CCSDT-unCP-atqz']['A24-4'], digits, 'CCSDT-unCP-atqz')
            qcdb.compare_values(-4.8102, h2kc * df['CCSDT-unCP-atz']['A24-4'], digits, 'CCSDT-unCP-atz')
            qcdb.compare_values(-4.4767, h2kc * df['MP2-CP-a56z']['A24-4'], digits, 'MP2-CP-a56z')
            qcdb.compare_values(-4.4282, h2kc * df['MP2-CP-a5z']['A24-4'], digits, 'MP2-CP-a5z')
            qcdb.compare_values(-4.4492, h2kc * df['MP2-CP-a6z']['A24-4'], digits, 'MP2-CP-a6z')
            qcdb.compare_values(-4.3078, h2kc * df['MP2-CP-adtz']['A24-4'], digits, 'MP2-CP-adtz')
            qcdb.compare_values(-3.9377, h2kc * df['MP2-CP-adz']['A24-4'], digits, 'MP2-CP-adz')
            qcdb.compare_values(-4.4721, h2kc * df['MP2-CP-aq5z']['A24-4'], digits, 'MP2-CP-aq5z')
            qcdb.compare_values(-4.3816, h2kc * df['MP2-CP-aqz']['A24-4'], digits, 'MP2-CP-aqz')
            qcdb.compare_values(-4.4646, h2kc * df['MP2-CP-atqz']['A24-4'], digits, 'MP2-CP-atqz')
            qcdb.compare_values(-4.1941, h2kc * df['MP2-CP-atz']['A24-4'], digits, 'MP2-CP-atz')
            qcdb.compare_values(-4.4833, h2kc * df['MP2-unCP-a56z']['A24-4'], digits, 'MP2-unCP-a56z')
            qcdb.compare_values(-4.5891, h2kc * df['MP2-unCP-a5z']['A24-4'], digits, 'MP2-unCP-a5z')
            qcdb.compare_values(-4.5411, h2kc * df['MP2-unCP-a6z']['A24-4'], digits, 'MP2-unCP-a6z')
            qcdb.compare_values(-4.7573, h2kc * df['MP2-unCP-adtz']['A24-4'], digits, 'MP2-unCP-adtz')
            qcdb.compare_values(-4.6150, h2kc * df['MP2-unCP-adz']['A24-4'], digits, 'MP2-unCP-adz')
            qcdb.compare_values(-4.5823, h2kc * df['MP2-unCP-aq5z']['A24-4'], digits, 'MP2-unCP-aq5z')
            qcdb.compare_values(-4.6374, h2kc * df['MP2-unCP-aqz']['A24-4'], digits, 'MP2-unCP-aqz')
            qcdb.compare_values(-4.6023, h2kc * df['MP2-unCP-atqz']['A24-4'], digits, 'MP2-unCP-atqz')
            qcdb.compare_values(-4.6749, h2kc * df['MP2-unCP-atz']['A24-4'], digits, 'MP2-unCP-atz')
            qcdb.compare_values(-3.8207, h2kc * df['HF-CP-a5z']['A24-4'], digits, 'HF-CP-a5z')
            qcdb.compare_values(-3.8217, h2kc * df['HF-CP-a6z']['A24-4'], digits, 'HF-CP-a6z')
            qcdb.compare_values(-3.7558, h2kc * df['HF-CP-adz']['A24-4'], digits, 'HF-CP-adz')
            qcdb.compare_values(-3.8160, h2kc * df['HF-CP-aqz']['A24-4'], digits, 'HF-CP-aqz')
            qcdb.compare_values(-3.7422, h2kc * df['HF-CP-atz']['A24-4'], digits, 'HF-CP-atz')
            qcdb.compare_values(-3.8286, h2kc * df['HF-unCP-a5z']['A24-4'], digits, 'HF-unCP-a5z')
            qcdb.compare_values(-3.8226, h2kc * df['HF-unCP-a6z']['A24-4'], digits, 'HF-unCP-a6z')
            qcdb.compare_values(-3.9955, h2kc * df['HF-unCP-adz']['A24-4'], digits, 'HF-unCP-adz')
            qcdb.compare_values(-3.8704, h2kc * df['HF-unCP-aqz']['A24-4'], digits, 'HF-unCP-aqz')
            qcdb.compare_values(-3.8598, h2kc * df['HF-unCP-atz']['A24-4'], digits, 'HF-unCP-atz')

            #qcdb.compare_values(-0.0394, h2kc * df['CCSDT-CP-adz']['A24-4'], digits, 'CCSDT-CP-adz')  # TODO this is delta_MP2^CCSD(T)
            #qcdb.compare_values(-0.1062, h2kc * df['CCSDT-CP-atz']['A24-4'], digits, 'CCSDT-CP-atz')  # TODO this is delta_MP2^CCSD(T)
            qcdb.compare_values(-4.5040, h2kc * df['CCSDT-CP-atqzadz']['A24-4'], digits, 'CCSDT-CP-atqzadz')
            qcdb.compare_values(-4.5708, h2kc * df['CCSDT-CP-atqzatz']['A24-4'], digits, 'CCSDT-CP-atqzatz')
            qcdb.compare_values(-4.5115, h2kc * df['CCSDT-CP-aq5zadz']['A24-4'], digits, 'CCSDT-CP-aq5zadz')
            qcdb.compare_values(-4.5783, h2kc * df['CCSDT-CP-aq5zatz']['A24-4'], digits, 'CCSDT-CP-aq5zatz')

            #qcdb.compare_values(-0.0915, h2kc * df['CCSDT-ave-adz']['A24-4'], digits, 'CCSDT-ave-adz')  # TODO this is delta_MP2^CCSD(T)
            #qcdb.compare_values(-0.1208, h2kc * df['CCSDT-ave-atz']['A24-4'], digits, 'CCSDT-ave-atz')  # TODO this is delta_MP2^CCSD(T)
            qcdb.compare_values(-4.5335, h2kc * df['MP2-ave-atqz']['A24-4'], digits, 'MP2-ave-atqz')
            qcdb.compare_values(-4.5272, h2kc * df['MP2-ave-aq5z']['A24-4'], digits, 'MP2-ave-aq5z')
            qcdb.compare_values(-4.6250, h2kc * df['CCSDT-ave-atqzadz']['A24-4'], digits, 'CCSDT-ave-atqzadz')
            qcdb.compare_values(-4.6542, h2kc * df['CCSDT-ave-atqzatz']['A24-4'], digits, 'CCSDT-ave-atqzatz')
            qcdb.compare_values(-4.6187, h2kc * df['CCSDT-ave-aq5zadz']['A24-4'], digits, 'CCSDT-ave-aq5zadz')
            qcdb.compare_values(-4.6480, h2kc * df['CCSDT-ave-aq5zatz']['A24-4'], digits, 'CCSDT-ave-aq5zatz')

    except KeyError, e:
        print e
        pass
    else:
        print 'End of Tests'
