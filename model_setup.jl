# Set up of the reaction network model using Catalyst.jl
using Catalyst 

growth_rtc_model = @reaction_network begin
    # Import nutrient                                   
    imp(et, vt, s, kt),                                                 ∅ => si

    # Metabolization nutrient
    met(em, vm, si, km),                                                si => ∅
    met(em, vm, si, km) * ns,                                           ∅ => e

    # Transcription
    tx(wri, e, thetari),                                                ∅ => mri
    tx(we, e, thetatx),                                                 ∅ => mt 
    tx(we, e, thetatx),                                                 ∅ => mm
    tq(wq, e, thetatx, q, kq, hq),                                      ∅ => mq
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

    # Antibiotic influx/eflux
    #(pin, pout),                                                        abx <--> abxi        
    pin * abx,                                                                ∅ --> abxi # before abx
    pout * abxi,                                                              abxi --> ∅ # before abx

    # Ribosome binding antibiotics

    #(kon, koff),                                                        abxi + cri <--> zmri
    kon,                                                                abxi + cri --> zmri
    koff,                                                               zmri --> abxi + cri

    #(kon, koff),                                                        abxi + ct <--> zmt
    kon,                                                                abxi + ct --> zmt
    koff,                                                               zmt --> abxi + ct

    #(kon, koff),                                                        abxi + cm <--> zmm
    kon,                                                                abxi + cm --> zmm
    koff,                                                               zmm --> abxi + cm

    #(kon, koff),                                                        abxi + cq <--> zmq
    kon,                                                                abxi + cq --> zmq
    koff,                                                               zmq --> abxi + cq
   
    #(kon, koff),                                                        abxi + cr <--> zmr
    kon,                                                                abxi + cr --> zmr
    koff,                                                               zmr --> abxi + cr

    #(kon, koff),                                                        abxi + ca <--> zma
    kon,                                                                abxi + ca --> zma
    koff,                                                               zma --> abxi + ca

    #(kon, koff),                                                        abxi + cb <--> zmb
    kon,                                                                abxi + cb --> zmb
    koff,                                                               zmb --> abxi + cb
 
    # Damage rate -> scale antibiotic concentration with a factor?
    abxi * kdam,                                                         rh --> rd
    # Damage of ribosome complex
    abxi * kdam,                                                         (cri, ct, cm, cq, cr, ca, cb) --> (crid, ctd, cmd, cqd, crd, cad, cbd)
    # Damage of ribosome-antibiotic complex 
    abxi * kdam,                                                         (zmri, zmt, zmm, zmq, zmr, zma, zmb) --> (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd)


    # Degradation of damaged RNA
    kdeg,                                                               rd --> ∅
    kdeg,                                                               (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd) --> ∅
    kdeg,                                                               (crid, ctd, cmd, cqd, crd, cad, cbd) --> ∅

    # Tagging
    vtag(a, rd, ktag, ka),                                              rd => rt
    vtag(a, rd, ktag, ka),                                              (zmrid, zmtd, zmmd, zmqd, zmrd, zmad, zmbd) => (zmrit, zmtt, zmmt, zmqt, zmrt, zmat, zmbt)
    vtag(a, rd, ktag, ka),                                              (crid, ctd, cmd, cqd, crd, cad, cbd) => (crit, ctt, cmt, cqt, crt, cat, cbt)

    # Repair
    vrep(b, rt, krep, kb),                                              rt => rh 
    vrep(b, rt, krep, kb),                                              (zmrit, zmtt, zmmt, zmqt, zmrt, zmat, zmbt) => (zmri, zmt, zmm, zmq, zmr, zma, zmb)
    vrep(b, rt, krep, kb),                                              (crit, ctt, cmt, cqt, crt, cat, cbt) => (cri, ct, cm, cq, cr, ca, cb )

    # Dilution add damage complexes!
    lam(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr, m),               (e, si, abxi, mri, mt, mm, mq, ma, mb, mr, et, em, q, a, b, r, rh, cri, ct, cm, cq, ca, cb, cr, zmri, zmt, zmm, zmq, zma, zmb, zmr, rd, crid, ctd, cmd, cqd, cad, cbd, crd, zmrid, zmtd, zmmd, zmqd, zmad, zmbd, zmrd, rt, crit, ctt, cmt, cqt, cat, cbt, crt, zmrit, zmtt, zmmt, zmqt, zmat, zmbt, zmrt) --> ∅
    
    # Energy consumption during translation
    econ(e, gmax, Kgamma, cri, ct, cm, cq, ca, cb, cr),                 e => ∅
end