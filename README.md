# ***TITLE OF PAPER***

## Description

In our previous work, we presented two different optimization schemes that estimate parameters for a Partial Differential Equation (PDE) model of cancer invasion based on synthetic spatial 1D data[^1]. Given the positive results with the spatial 1D dataset on variable densities, we now extend the study to spatial 2D datasets, which requires discretizing the PDE model in one more dimension within its numerical solver. We now focus on a spatial 2D individual-based model (IBM) that describes the micro-scale movements and allows us to study migration and invasion at the level of individual cells. Within this new modelling framework, the simulated patterns can be compared with cancer invasion patterns observed in in vitro or ex vivo assays, such as organotypic cultures, enabling us to draw parameter estimates on a more realistic basis. 

This repository contains all the simulation code and results of parameter optimizations on two different sets of invasion patterns recorded in organotypic cultures: the first describes the invasion of Squamous Cells Carcinoma (SCC) in an in vitro organotypic culture (mainly contains synthetic stroma composed of collagen gel embedded with fibroblasts) over 14 days [^2], whilst the second demonstrates the invasion of a T98G glioma spheroid in an ex vivo rat brain slices culture over a period of 3 days [^3].

[^1]: Xiao, Y., Thomas, L. & Chaplain, M. (2021), ‘Calibrating models of cancer invasion and metastasis: parameter inference using approximate Bayesian computation and gradient matching.’, R. Soc. Open Sci. 8(202237).
[^2]: Nystrom, S., Thomas, G., Stone, M., Mackenzie, I., Hart, I. & Marshall, J. (2005), ‘Development of a quantitative method to analyse tumour cell invasion in organotypic culture’, J. Pathol 205, 468–475.
[^3]: Matsumura, H., Ohnishi, T., Kanemura, Y., Maruno, M. & Yoshimine, T. (2000), ‘Quantitative analysis of glioma cell invasion by confocal laser scanning microscopy in a novel brain slice model’, Biochem. Biophys. Res. Commun. 269, 513–520.
