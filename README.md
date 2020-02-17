# SpeckledNeurons_control

---

**a lot of code need to be cleaned and properly commented**

---

main files:
* **calibrate_beads_sample.m** - calibrate the bottom camera with the DMD, and create a T.mat file which must be called afterwards since gives the parameters of the affine transformation between the two spaces
* **thaFancyGameOfSearchingBeads.m** - search the beads and helps to setup the recordings
* **GenerateAndRecordCaActBeads.m**  - main script, which takes beads centroids collected with "thaFancyGameOfSearchingBeads" and do the needed stuff to excite them with the temporal activity using simple exponentials, or taking some random piece of electrophisiologically recorded data

External projects used within this code:

* **DMDConnect-master**
 https://github.com/deichrenner/DMDConnect

* **MLspike pieces**
https://github.com/MLspike
from the paper "Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo." Deneux et al, Nat. Comm.  
  **brick-master/**  
  **spikes-master/**  

* **spikefinder.train**
datasets from http://spikefinder.codeneuro.org/
