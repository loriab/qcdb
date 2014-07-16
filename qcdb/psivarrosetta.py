
useme2psivar = {
    #'usemeraw': 'HF TOTAL ENERGY',
    'usemeraw': 'SCF TOTAL ENERGY',

    'mp2.usemecorl': 'MP2 CORRELATION ENERGY',
    'mp3.usemecorl': 'MP3 CORRELATION ENERGY',
    'mp4.usemecorl': 'MP4 CORRELATION ENERGY',
    'ccsd.usemecorl': 'CCSD CORRELATION ENERGY',
    'ccsdt.usemecorl': 'CCSD(T) CORRELATION ENERGY',
    'ccsdfullt.usemecorl': 'CCSDT CORRELATION ENERGY',
    'ccsdtq.usemecorl': 'CCSDT(Q) CORRELATION ENERGY',
    
    'ccsdt.usemecrct': '(T) CORRECTION ENERGY',
    'ccsdtq.usemecrct': '(Q) CORRECTION ENERGY',

    'mp2.usemetrip': 'MP2 SAME-SPIN CORRELATION ENERGY',

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

