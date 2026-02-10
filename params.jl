# Definition of the parameters and initial conditions for the model 
sf = (1e6 * 1 /(6.02e23 * 1e-15)) 

kbi_value =  1 / sf # [cell /min molecs], Growth
#Kgamma_growth = 7 * sf # [molecs /cell], Growth
km_value = 1000 * sf # [molecs /cell], Growth
kon_value = 0.023 
kon_value_test = 2
kq_value = 152219 * sf # [molecs /cell], Growth
kt_value = 1000 #* sf # [molecs], Growth 
m_value = 1e8 * sf # aa, Growth
thetari_value = 426.87 * sf	# [molecs /cell], Growth 
thetatx_value = 4.38 * sf # [molecs /cell], Growth 
we_value = 4.14 * sf # molecs /min cell, Growth
wri_value = 930 * sf # molecs /min cell, Growth
wq_value = 948.93 * sf #molecs /min cell
s_value = 1e4 #* sf


params = Dict(
    :abx => 0, # µM -> Change
    :c => 0.01, # Unitless, Rtc
    :dm => 0.1, # min−1, Growth
    :gmax => 1260, # aa /min molecs
    :hq => 4, # Unitless, Growth
    :ka => 20, # µM, Rtc
    :kb => 16, # µM, Rtc
    :kbi => kbi_value, # (min µM)−1, Growth
    #:kcm => 0.00599, # (min µM)−1, Growth
    :kdeg => 1, # min−1, Rtc
    :kdam => 0,#1, #?
    :kdiss => 0.006, # min−1, Rtc
    :Kgamma => 255.73, # µM, 
    :kh => 250, # µM, Rtc
    :km => km_value, # µM, Growth
    :koff => 2.76 * kon_value, # KD in [min-1]
    :kon => kon_value, # µM-1 min-1, Antibiotics
    :kq => kq_value, # µM, Growth
    :kr => 0.125, # µM, Rtc
    :krep => 0,#15.67, # min−1, Rtc
    :kt => kt_value, # µM, Growth
    :ktag => 0,#9780, # min−1, Rtc
    :kub => 1, # min−1, Growth
    :l => 50, # Unitless, Rtc
    :m => m_value, # aa, Growth
    :na => 338, # aa, Rtc
    :nb => 408, # aa, Rtc
    :nr => 532 * 6, # aa, Rtc
    :nri => 7459, # aa/molecs, Growth
    :ns => 0.5, # Unitless, Growth
    :nx => 300, # aa/molecs, Growth 
    :pin => 1.48, # min-1, Antibiotics
    :pout => 0.4, # min-1, Antibiotics
    :s => s_value, # µM, Growth
    :thetari => thetari_value,  # µM, Growth
    :thetatx => thetatx_value, # µM, Growth
    :vm	=> 5800, # min−1, Growth 
    :vmax => 39.51, # min−1, Rtc
    :vt => 726, # min−1, Growth
    :wba => 0,#5e−5, # min−1, Rtc
    :we => we_value, # µM·min−1, Growth
    :wr => 0,#1e−6, # µM·min−1, Rtc
    :wri => wri_value, # µM·min−1, Growth
    :wq	=> wq_value, # µM·min−1, Growth
)

e_value = 1000.0 * sf # molecs /cell?, Growth
rh_value = 10 * sf
rh_value_r = 11.29

u0 = Dict(
    :ma => 0.0, :mb => 0.0, :mm => 0.0, :mq => 0.0, :mr => 0.0, :mri => 0.0, :mt => 0.0, # mRNAs
    :a => 0.0, :b => 0.0, :em => 0.0, :et => 0.0, :q => 0.0, :r => 0.0, :rh => rh_value, # Proteins, rh (growth)
    :ca => 0.0, :cb => 0.0, :cm => 0.0, :cq => 0.0, :cr => 0.0, :cri => 0.0, :ct => 0.0, # Ribosome-mRNA complexes
    :cad => 0.0, :cbd => 0.0, :cmd => 0.0, :cqd => 0.0, :crd => 0.0, :crid => 0.0, :ctd => 0.0, # Damaged ribosome-mRNA complexes
    :cat => 0.0, :cbt => 0.0, :cmt => 0.0, :cqt => 0.0, :crt => 0.0, :crit => 0.0, :ctt => 0.0, # Tagged ribosome-mRNA complexes
    :rd => 0.0, :rt => 0.0, # Damaged ribosome
    :zma => 0.0, :zmb => 0.0, :zmm => 0.0, :zmq => 0.0, :zmr => 0.0, :zmri => 0.0, :zmt => 0.0, # Antibiotic-ribosome-mRNA (zombie) complex
    :zmad => 0.0, :zmbd => 0.0, :zmmd => 0.0, :zmqd => 0.0, :zmrd => 0.0, :zmrid => 0.0, :zmtd => 0.0, # Damaged zombie complex
    :zmat => 0.0, :zmbt => 0.0, :zmmt => 0.0, :zmqt => 0.0, :zmrt => 0.0, :zmrit => 0.0, :zmtt => 0.0, # Tagged zombie complex
    :abxi => 0.0, # µM
    :e => e_value, # µM
    :si => 0.0 # µM
)