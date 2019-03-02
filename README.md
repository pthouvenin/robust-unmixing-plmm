# A hierarchical Bayesian model accounting for endmember variability and abrupt spectral changes to unmix multitemporal hyperspectral images

**Description:** Matlab codes associated with the method described in 

>P.-A. Thouvenin, N. Dobigeon and J.-Y. Tourneret - <strong>A hierarchical Bayesian model accounting for endmember variability and abrupt spectral changes to unmix multitemporal hyperspectral images</strong>, <em>IEEE Trans. Comput. Imag.</em>, vol. 4, no. 1, pp. 32-45, Mar. 2018.

**Author:** P.-A. Thouvenin, pierreantoine[dot]thouvenin[at]gmail[dot]com

**Experiments:** to run a representative example of the real data experiments reported in the article, configure and run the `main_real_data.m` script. The script `main_extract_data.m` can be used to extract the hyperspectral data from the raw data file included in the `data/raw_data` folder. A `.mat` file obtained after data extraction is already provided in the `data` folder.

**Dependencies:** the present codes includes MATLAB functions described in the following publications, and developed by their respective authors.

> [1] J. M. Nascimento and J. M. Bioucas-Dias - <strong>Vertex component analysis: a fast algorithm to unmix hyperspectral data</strong>, <em>IEEE Trans. Geosci. Remote Sens.</em>, vol. 43, no. 4, pp. 898--910, Apr. 2005.

> [2] J. M. Bioucas-Dias and M. A. T. Figueiredo - <strong>Alternating direction algorithms for constrained sparse regression: Application to hyperspectral unmixing</strong>, <em>Proc. IEEE GRSS Workshop Hyperspectral Image Signal Process.: Evolution in Remote Sens. (WHISPERS).</em>, Reykjavik, Iceland, Jun. 2010.

> [3] J. Bioucas-Dias and J. Nascimento - <strong>Hyperspectral subspace identification</strong>, <em>IEEE Transactions on Geoscience and Remote Sensing.</em>, vol. 46, no. 8, pp. 2435-2445, 2008.

> [4] J. M. Bioucas-Dias - <strong>A variable splitting augmented Lagrangian approach to linear spectral unmixing</strong>, <em>Proc. IEEE GRSS Workshop Hyperspectral Image Signal Process.: Evolution in Remote Sens. (WHISPERS).</em>, Grenoble, France, Aug. 2009.

> [5] V. Mazet, - <strong>Simulation d'une distribution gaussienne tronquée sur un intervalle fini</strong>, Technical Report, Université de Strasbourg/CNRS, 2012. [<a href="http://miv.u-strasbg.fr/mazet/rtnorm/">Code</a>]

> [6] J. Tursa, MTIMESX - Fast Matrix Multiply with Multi-Dimensional Support, [<a href="https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support">Code on Matlab FileExchange</a>], [<a href="https://github.com/cybertk/mtimesx">Code on Github</a>]
