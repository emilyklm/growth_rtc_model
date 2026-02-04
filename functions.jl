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

