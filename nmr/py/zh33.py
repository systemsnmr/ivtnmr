# Zoom UP1-H33 in 2D 15N-HSQC 

#7.422 / 116.310
#6.690 / 116.807
#7.351 / 118.887

# R82, L16, T25
# 7.582 / 120.642
# 9.628 / 126.098
# 7.351 / 118.887

H = 7.422
N = 116.310
H_min = H-0.1
H_max = H+0.1
N_min = N-1.2
N_max = N+1.2

fullrange = putil.DataChecks.getNMRDataOfSelectedFramePrintMsg().getFullPhysicalRange()
newRange = fullrange
newRange[0].setStart(float( H_max ))
newRange[0].setEnd(float( H_min ))
newRange[0].setUnit("ppm")
newRange[1].setStart(float( N_max ))
newRange[1].setEnd(float( N_min ))
newRange[1].setUnit("ppm")
XCMD(".zx", 1, newRange)