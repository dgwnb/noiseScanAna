make ./bin/analyzeNoiseScanLBLQ1
# ./bin/analyzeNoiseScanLBLQ1 trigRtTrigCnt Input/scanList/trigRtTrigCnt
# ./bin/analyzeNoiseScanLBLQ1 reverse Input/scanList/reverse

./bin/analyzeNoiseScanLBLQ1 USBpixXC Input/scanList/USBpixXC

Noise scan on FE1/LBLQ1 is done using USBpix/STcontrol, which can be compared to the results from RCE/calibGui.

From the noise occupancy distributions, one can see in general the two systems see very consistent picture: At lowest threshold scanned (AltFine=60) noise burst happened, which is then removed with increasing threshold, but the right half of the chip is still more noisy than the left half. With the the threshold increased further the left-right asymmetry is also removed, and the noise in the ganged pixel region also gradually disappears.

The difference in the noise occupancy seen by two systems is that USBpix observes more noisy ganged pixel region than RCE at low threshold.

OccupancyAnimation_FE1USBpix.gifOccupancyAnimation_FE1RCE.gif

To make fair comparison between two set of results and focus on the random noise, the following mask is used to remove ganged pixels and hot pixels. NB it is different from the masks have been using to analyze noise scan from RCE so far, which remove pixels keep firing in 10 adjacent scans (thus the results are not comparable)

MaskFE1.gif

The noise occupancy changing as a function of threshold is summarized below. It is noted that though the two curves have similar behavior

