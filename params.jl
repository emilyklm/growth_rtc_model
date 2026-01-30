# Definition of the parameters and initial conditions for the model 
sf = (1e6 * 1 /(6.02e23 * 1e-15))

# vtag_val = 3254 * sf # uM 
kbi_value =  1 * sf # [cell /min molecs], Growth
# Kgamma = 7 * sf # [molecs /cell], Growth
km_value = 1000 * sf # [molecs /cell], Growth
kq_value = 152219 * sf # [molecs /cell], Growth
kt_value = 1000 * sf # [molecs], Growth
s_value = 1e4 * sf # molecules, Growth 
thetari_value = 426.87 * sf	# [molecs /cell], Growth 
thetatx_value = 4.38 * sf # [molecs /cell], Growth 



params = Dict(
# :vtag => vtag_val
:c => 0.01, # Unitless, Rtc
:dm => 0.1, # min−1, Growth
:hq => 4, # Unitless, Rtc
:ka => 20, # µM, Rtc
:kb => 16, # µM, Rtc
:kbi => kbi_value, # (min µM)−1, Growth
:kcm => 0.00599, # (min µM)−1, Growth
:kdeg => 1, # min−1, Rtc
# :kdam ?
:kdiss => 0.006, # min−1, Rtc
:Kgamma => 255.73, # µM, 
:kh => 250, # µM, Rtc
:km => km_value, # µM, Growth
:kq => kq_value, # µM, Growth
:kr => 0.125, # µM, Rtc
:krep => 15.67, # min−1, Rtc
:kt => kt_value, # µM, Growth
:ktag => 9780, # min−1, Rtc
:kub => 1, # min−1, Growth
:l => 50, # Unitless, Rtc
:m => 1e8, # aa, Growth
:na => 338, # aa, Rtc
:nb => 408, # aa, Rtc
:nr => 532, # aa, Rtc
:nri => 7459, # aa/molecs, Growth
:ns => 0.5, # Unitless, Growth
:nx => 300, # aa/molecs, Growth 
:s => s_value, # µM, Growth
:thetari => thetari_value,  # µM, Growth
:thetatx => thetatx_value, # µM, Growth
:vm	=> 5800, # min−1, Growth 
:vmax => 39.51, # min−1, Rtc
:vt => 726, # min−1, Growth
:wba => 5·10−5, # min−1, Rtc
:we => 0.0069, # µM·min−1, Growth
:wr => 1·10−6, # µM·min−1, Rtc
:wri => 1.54, # µM·min−1, Growth
:wq	=> 1.58, # µM·min−1, Growth

u0 = Dict(

)