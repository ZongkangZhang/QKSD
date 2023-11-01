# QKSD
These codes are accompanied to the following paper:

  * Zongkang Zhang, Anbang Wang, Xiaosi Xu, Ying Li, Measurement-efficient quantum Krylov subspace diagonalisation, [arXiv:2301.13353](https://arxiv.org/abs/2301.13353).

## Requirement
These codes are written in Mathematica and have been tested in Mathematica 13. 

## Description
  * In the folder [measurement_overhead_benchmarking](https://github.com/ZongkangZhang/QKSD/tree/main/measurement_overhead_benchmarking): [chain&ladder.nb](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/chain%26ladder.nb) and [random_graph.nb](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/random_graph.nb) calculate the measurement overhead approaching the accuracy of classical Lanczos algorithm for different quantum Krylov subspace diagonalisation (KSD) algorithms. The corresponding results are collected in the subfolder [data](https://github.com/ZongkangZhang/QKSD/tree/main/measurement_overhead_benchmarking/data), which is used in [empirical_distribution.nb](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/empirical_distribution.nb) to plot the empirical distribution of the measurement overhead. Taking the instance (Heisenberg, chain, d=5) as an example, [example.nb](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/example.nb) gives a comparison between different quantum KSD algorithms. In the above scripts, packages [Qubits_package.m](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/Qubits_package.m) and [ExactKrylov_package.m](https://github.com/ZongkangZhang/QKSD/blob/main/measurement_overhead_benchmarking/ExactKrylov_package.m) are imported.
  * In the folder [practical_implementation](https://github.com/ZongkangZhang/QKSD/tree/main/necessary_measurement_number): [regularisation.nb](https://github.com/ZongkangZhang/QKSD/blob/main/necessary_measurement_number/regularisation.nb) shows that the sufficient measurement costs computed by the theorem are close to the necessary measurement costs when using the regularisation method. [thresholding.nb](https://github.com/ZongkangZhang/QKSD/blob/main/necessary_measurement_number/thresholding.nb) compares the necessary measurement costs between different quantum KSD algorithms when using thresholding procedure. [Gaussian-power.nb](https://github.com/ZongkangZhang/QKSD/blob/main/necessary_measurement_number/Gaussian-power.nb) illustrates the Gaussian-power basis is a filter, which effects the choice of tau. The above scripts import the package [QLanczos_package.m](https://github.com/ZongkangZhang/QKSD/blob/main/necessary_measurement_number/QLanczos_package.m) and take the instance (Heisenberg, chain, d=5) as an example.
  * [ratio.nb](https://github.com/ZongkangZhang/QKSD/blob/main/ratio.nb) computes the ratio between cost upper bound and norm upper bound of the Gaussian-power basis.

## License
This code is under the [MIT license](https://opensource.org/licenses/MIT).
