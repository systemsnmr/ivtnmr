# Zoom imino region of SMN

fullrange = putil.DataChecks.getNMRDataOfSelectedFramePrintMsg().getFullPhysicalRange()
newRange = fullrange
newRange[0].setStart(float( 14.8 ))
newRange[0].setEnd(float( 10.3 ))
newRange[0].setUnit("ppm")
XCMD(".zx", 1, newRange)