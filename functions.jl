# Defined functions

# Nutrient import function
imp(et, vt, s, kt) = et * vt * s / (kt + s)

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
function txab(vmax, kh, kdiss, r, l, c, rt, kr, wx, e, thetatx)         # e is not a parameter anymore
    voc = (vmax * e) / (kh + e)                                           # Conversion rate of sigma
    sigmao = (active_rtcr(r, l, c, rt, kr) * voc) / kdiss
    return sigmao * tx(wx, e, thetatx)
end

# Gamma function
gamma(e, gmax, Kgamma) = (gmax * e) / (Kgamma + e)

# Translation rate function for all translations -> Add rtc mRNA to buid a ribosome-RNA complex! 
tl(cx, e, gmax, Kgamma, nx) = cx * (gamma(e, gmax, Kgamma) / nx)

# Energy consumption during translation
econ(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr) = gamma(e, gmax, Kgamma) * (cri + ct + cm + cq + ca + cb + cr) 

# Tagging rate 
vtag(a, rd, ktag, ka) = (a * rd * ktag) / (rd + ka)

# Repair rate
vrep(b, rt, krep, kb) = (b * rt * krep) / (rt + kb)

# Growth function
lam(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr, m) = (gamma(e, gmax, Kgamma) * (cri + ct + cm + cq + ca + cb + cr)) / m

# Antibiotic binding
bin(cx, abxi, kon) = cx * abxi * kon

# Antibiotic unbinding
unbin(zmx, koff) = zmx * koff

# Other functions
function with_drug(params, drug_params, drug::Symbol; abx=0.0)
    p = deepcopy(params)
    dp = drug_params[drug]
    p[:abx]  = abx
    p[:pin]  = dp[:pin]
    p[:pout] = dp[:pout]
    p[:kon]  = dp[:kon]
    p[:koff] = dp[:koff]

    return p
end
