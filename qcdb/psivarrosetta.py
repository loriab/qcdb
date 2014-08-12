
useme2psivar = {

    # <<<  DFT  >>>

    'b3lyp.usemeraw': 'B3LYP FUNCTIONAL TOTAL ENERGY',
    'b3lypd2.usemedash': 'B3LYP-D2 DISPERSION CORRECTION ENERGY',
    'b3lypd3.usemedash': 'B3LYP-D3 DISPERSION CORRECTION ENERGY',
    'b3lypd3bj.usemedash': 'B3LYP-D3BJ DISPERSION CORRECTION ENERGY',
    'b3lypxdm.usemedash': 'B3LYP-XDM DISPERSION CORRECTION ENERGY',

    'b2plyp.usemeraw': 'B2PLYP TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def
    'b2plypd2.usemedash': 'B2PLYP-D2 DISPERSION CORRECTION ENERGY',
    'b2plypd3.usemedash': 'B2PLYP-D3 DISPERSION CORRECTION ENERGY',
    'b2plypd3bj.usemedash': 'B2PLYP-D3BJ DISPERSION CORRECTION ENERGY',

    'b970.usemeraw': 'B970 FUNCTIONAL TOTAL ENERGY',
    'b970d2.usemedash': 'B970-D2 DISPERSION CORRECTION ENERGY',

    'b97.usemeraw': 'B97 FUNCTIONAL TOTAL ENERGY',
    'b97d2.usemedash': 'B97-D2 DISPERSION CORRECTION ENERGY',
    'b97d3.usemedash': 'B97-D3 DISPERSION CORRECTION ENERGY',
    'b97d3bj.usemedash': 'B97-D3BJ DISPERSION CORRECTION ENERGY',

    'bp86.usemeraw': 'BP86 FUNCTIONAL TOTAL ENERGY',
    'bp86d2.usemedash': 'BP86-D2 DISPERSION CORRECTION ENERGY',
    'bp86d3.usemedash': 'BP86-D3 DISPERSION CORRECTION ENERGY',
    'bp86d3bj.usemedash': 'BP86-D3BJ DISPERSION CORRECTION ENERGY',

    'wb97x.usemeraw': 'WB97X FUNCTIONAL TOTAL ENERGY',
    'wb97xd.usemeraw': 'WB97X-D TOTAL ENERGY',
    'wb97xd.usemedash': 'WB97X-D DISPERSION CORRECTION ENERGY',

    'wb97x2.usemeraw': 'WB97X-2 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'm052x.usemeraw': 'M05-2X FUNCTIONAL TOTAL ENERGY',
    'm052xd3.usemedash': 'M05-2X-D3 DISPERSION CORRECTION ENERGY',

    'm062x.usemeraw': 'M06-2X FUNCTIONAL TOTAL ENERGY',
    'm062xd3.usemedash': 'M06-2X-D3 DISPERSION CORRECTION ENERGY',

    'pbe.usemeraw': 'PBE FUNCTIONAL TOTAL ENERGY',
    'pbed2.usemedash': 'PBE-D2 DISPERSION CORRECTION ENERGY',
    'pbed3.usemedash': 'PBE-D3 DISPERSION CORRECTION ENERGY',
    'pbed3bj.usemedash': 'PBE-D3BJ DISPERSION CORRECTION ENERGY',

    'pbe0.usemeraw': 'PBE0 FUNCTIONAL TOTAL ENERGY',
    'pbe0d2.usemedash': 'PBE0-D2 DISPERSION CORRECTION ENERGY',
    'pbe0d3.usemedash': 'PBE0-D3 DISPERSION CORRECTION ENERGY',
    'pbe0d3bj.usemedash': 'PBE0-D3BJ DISPERSION CORRECTION ENERGY',

    'xyg3.usemeraw': 'XYG3 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'vv10.usemeraw': 'VV10 FUNCTIONAL TOTAL ENERGY',

    'lcvv10.usemeraw': 'LC-VV10 FUNCTIONAL TOTAL ENERGY',

    'dsdpbep86.usemeraw': 'DSD-PBEP86 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'm08hx.usemeraw': 'M08-HX FUNCTIONAL TOTAL ENERGY',
    'm08so.usemeraw': 'M08-SO FUNCTIONAL TOTAL ENERGY',
    'm11.usemeraw': 'M11 FUNCTIONAL TOTAL ENERGY',
    'm11l.usemeraw': 'M11L FUNCTIONAL TOTAL ENERGY',

    'pbe02.usemeraw': 'PBE0-2 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'dldf.usemeraw': 'DLDF FUNCTIONAL TOTAL ENERGY',
    'dldfd.usemedash': 'DLDF+D DISPERSION CORRECTION ENERGY',

    # <<<  WFN  >>>

    #'usemeraw': 'HF TOTAL ENERGY',
    'usemeraw': 'SCF TOTAL ENERGY',

    'mp2.usemecorl': 'MP2 CORRELATION ENERGY',
    'mp3.usemecorl': 'MP3 CORRELATION ENERGY',
    'mp4.usemecorl': 'MP4 CORRELATION ENERGY',
    'ccsd.usemecorl': 'CCSD CORRELATION ENERGY',
    'ccsdt.usemecorl': 'CCSD(T) CORRELATION ENERGY',
    'ccsdfullt.usemecorl': 'CCSDT CORRELATION ENERGY',
    'ccsdtq.usemecorl': 'CCSDT(Q) CORRELATION ENERGY',

    'fno.usemecrct': 'FNO CORRECTION ENERGY',
    'fnoccsdt.usemecorl': 'CCSD(T) FNO CORRELATION ENERGY',

    'ccsdt.usemecrct': '(T) CORRECTION ENERGY',
    'ccsdtq.usemecrct': '(Q) CORRECTION ENERGY',

    'mp2.usemetrip': 'MP2 SAME-SPIN CORRELATION ENERGY',
    'ccsd.usemetrip': 'CCSD SAME-SPIN CORRELATION ENERGY',

    # <<<  F12  >>>

    'f12.usemeraw': 'HF-CABS TOTAL ENERGY',
    'mp2f12.usemecorl': 'MP2-F12 CORRELATION ENERGY',
    'ccsdaf12.usemecorl': 'CCSD-F12A CORRELATION ENERGY',
    'ccsdbf12.usemecorl': 'CCSD-F12B CORRELATION ENERGY',
    'ccsdcf12.usemecorl': 'CCSD-F12C CORRELATION ENERGY',
    'ccsdnstaf12.usemecorl': 'CCSD(T)-F12A CORRELATION ENERGY',
    'ccsdstaf12.usemecorl': 'CCSD(T*)-F12A CORRELATION ENERGY',
    'ccsdtaf12.usemecorl': 'CCSD(T**)-F12A CORRELATION ENERGY',
    'ccsdnstbf12.usemecorl': 'CCSD(T)-F12B CORRELATION ENERGY',
    'ccsdstbf12.usemecorl': 'CCSD(T*)-F12B CORRELATION ENERGY',
    'ccsdtbf12.usemecorl': 'CCSD(T**)-F12B CORRELATION ENERGY',
    'ccsdnstcf12.usemecorl': 'CCSD(T)-F12C CORRELATION ENERGY',
    'ccsdstcf12.usemecorl': 'CCSD(T*)-F12C CORRELATION ENERGY',
    'ccsdtcf12.usemecorl': 'CCSD(T**)-F12C CORRELATION ENERGY',

    'ccsdnstabf12.usemecrct': '(T)-F12AB CORRECTION ENERGY',
    'ccsdstabf12.usemecrct': '(T*)-F12AB CORRECTION ENERGY',
    'ccsdtabf12.usemecrct': '(T**)-F12AB CORRECTION ENERGY',
    'ccsdnstcf12.usemecrct': '(T)-F12C CORRECTION ENERGY',
    'ccsdstcf12.usemecrct': '(T*)-F12C CORRECTION ENERGY',
    'ccsdtcf12.usemecrct': '(T**)-F12C CORRECTION ENERGY',

    'mp2f12.usemetrip': 'MP2-F12 SAME-SPIN CORRELATION ENERGY',
    'ccsdaf12.usemetrip': 'CCSD-F12A SAME-SPIN CORRELATION ENERGY',
    'ccsdbf12.usemetrip': 'CCSD-F12B SAME-SPIN CORRELATION ENERGY',
    'ccsdcf12.usemetrip': 'CCSD-F12C SAME-SPIN CORRELATION ENERGY',
    }

psivar2useme = dict((v, k) for k, v in useme2psivar.iteritems())

