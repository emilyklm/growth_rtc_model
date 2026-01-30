# Defined functions
# Nutrient import function
imp(et, vt, s, Kt) = et * vt * s / (Kt + s)

# Nutrient metabolism function
met(em, vm, si, Km) = em * vm * si / (Km + si)

# Transcription function for mr, mt, mm and rtcr -> wx is different for rtcr!, thetatx different for mr
tx(wx, e, thetatx) = wx * (e / (thetatx + e))

# Transcription function for mq 
tq(wx, e, thetatx, q, Kq, hq) = tx(wx, e, thetatx) * (1 / (1 + (q / Kq) ^hq))

# Active RtcR
function active_rtcr(r, l, c, rt, kr)
    alpha = rt / kr                                                              # Normalised concentration of rt 
    fa = (1 + alpha)^6 / (l * (1 + c * alpha)^6 + (1 + alpha)^6)                         # Fraction of open RtcR
    return fa * r
end

# Transcription function for rtca, rtcb, thetatx different for mr
function txab(vmax, kh, kdiss, r, l, c, rt, kr, wx, e, thetatx)         # a is not a parameter anymore
    voc = (vmax * e) / (kh + e)                                           # Conversion rate of sigma
    sigmao = (active_rtcr(r, l, c, rt, kr) * voc) / kdiss
    return sigmao * tx(wx, e, thetatx)
end

# Gamma function
gamma(e, gmax, Kgamma) = (gmax * e) / (Kgamma + e)

# Translation rate function for all translations -> Add rtc mRNA to buid a ribosome-RNA complex! kc * rh * mx -> cx from tl(nx, kc, rh, mx, gmax, atp, thetatl) = (1 / nx) * kc * rh * mx * ((gmax * atp) / (thetatl + atp)) 
tl(cx, e, gmax, Kgamma, nx) = cx * (gamma(e, gmax, Kgamma) / nx)

# Energy consumption during translation
econ(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr) = gamma(e, gmax, Kgamma) * (cri + ct + cm + cq + ca + cb + cr) 

# Tagging rate 
vtag(a, rd, ktag, ka) = (a * rd * ktag) / (rd + ka)

# Repair rate
vrep(b, rt, krep, kb) = (b * rt * krep) / (rt + kb)

# Growth function
lam(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr, m) = (gamma(e, gmax, Kgamma) * (cri + ct + cm + cq + ca + cb + cr)) / m


# Set up of the reaction network model using Catalyst.jl
using Catalyst 

growth_rtc_model = @reaction_network begin
    # Import nutrient                                   
    imp(et, vt, s, Kt),                                                 ∅ => si

    # Metabolization nutrient
    met(em, vm, si, Km),                                                si => ∅
    met(em, vm, si, Km) * ns,                                           ∅ => e

    # Transcription
    tx(wri, e, thetari),                                                ∅ => mri
    tx(we, e, thetatx),                                                 ∅ => mt 
    tx(we, e, thetatx),                                                 ∅ => mm
    tq(wq, e, thetatx, q, Kq, hq),                                      ∅ => mq
    tx(wr, e, thetatx),                                                 ∅ => mr
    txab(vmax, kh, kdiss, r, l, c, rt, kr, wba, e, thetatx),            ∅ => (ma, mb)

    # degradation -> or use only one?
    dm,                                                                (mri, mt, mm, mq, mr, ma, mb) --> ∅ 

    # ribosome binding of mRNA
    kbi,                                                                rh + mri --> cri
    kub,                                                                cri --> rh + mri

    kbi,                                                                rh + mt --> ct
    kub,                                                                ct --> rh + mt

    kbi,                                                                rh + mm --> cm
    kub,                                                                cm --> rh + mm

    kbi,                                                                rh + mq --> cq
    kub,                                                                cq --> rh + mq

    kbi,                                                                rh + mr --> cr
    kub,                                                                cr --> rh + mr

    kbi,                                                                rh + ma --> ca
    kub,                                                                ca --> rh + ma

    kbi,                                                                rh + mb --> cb
    kub,                                                                cb --> rh + mb

    # Translation 
    tl(cri, e, gmax, Kgamma, nri),                                      cri => rh + mri + rh
    tl(ct, e, gmax, Kgamma, nx),                                        ct => rh + mt + et
    tl(cm, e, gmax, Kgamma, nx),                                        cm => rh + mm + em
    tl(cq, e, gmax, Kgamma, nx),                                        cq => rh + mq + q
    tl(cr, e, gmax, Kgamma, nr),                                        cr => rh + mr + r
    tl(ca, e, gmax, Kgamma, na),                                        ca => rh + ma + a
    tl(cb, e, gmax, Kgamma, nb),                                        cb => rh + mb + b
    
    # no kin -> rh is created via translation

    # Ribosome binding antibiotics
    abx * kcm,                                                          cri --> zmri
    abx * kcm,                                                          ct --> zmt
    abx * kcm,                                                          cm --> zmm
    abx * kcm,                                                          cq --> zmq
    abx * kcm,                                                          cr --> zmr
    abx * kcm,                                                          ca --> zma
    abx * kcm,                                                          cb --> zmb
    
 
    # Damage rate -> scale antibiotic concentration with a factor?
    abx * kdam,                                                         rh --> rd
    # Or damage of ribosome complex
    abx * kdam,                                                         (cri, ct, cm, cq, cr, ca, cb) --> (crid, ctd, cmd, cqd, crd, cad, cbd)
    # Or damage of ribosme-antibiotic complex -> back to zombie complex
    abx * kdam,                                                         (zmri, zmt, zmm, zmq, zmr, zma, zmb) --> (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd)
    # Or damage of ribosme-antibiotic complex and release of antibiotic
    abx * kdam,                                                         (zmri, zmt, zmm, zmq, zmr, zma, zmb) --> (crid, ctd, cmd, cqd, crd, cad, cbd)

    # Degradation of damaged RNA
    kdeg,                                                               rd --> ∅
    kdeg,                                                               (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd) --> ∅
    kdeg,                                                               (crid, ctd, cmd, cqd, crd, cad, cbd) --> ∅

    # Tagging
    vtag(a, rd, ktag, ka),                                              rd => rt
    vtag(a, rd, ktag, ka),                                              (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd) => (zmrit, zmtt, zmmt, zmqt, zmrt, zmat, zmbt)
    vtag(a, rd, ktag, ka),                                              (crid, ctd, cmd, cqd, crd, cad, cbd) --> (crit, ctt, cmd, cqt, crt, cat, cbt)

    # Repair
    vrep(b, rt, krep, kb),                                              rt => rh 
    vrep(b, rt, krep, kb),                                              (zmrit, zmtt, zmmt, zmqt, zmrt, zmat, zmbt) => (zmri, zmt, zmm, zmq, zmr, zma, zmb)
    vrep(b, rt, krep, kb),                                              (crit, ctt, cmd, cqt, crt, cat, cbt) => (cri, ct, cm, cq, cr, ca, cb )

    # Dilution add damage complexes!
    lam(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr, m),               (rh, et, em, q, a, b, r, si, mri, mt, mm, mq, ma, mb, mr, e, cri, ct, cm, cq, cs, cb, cr, rd, rt, zmri, zmt, zmm, zmq, zma, zmb, zmr) --> ∅
    
    # Energy consumption during translation
    econ(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr)                  e => ∅
end