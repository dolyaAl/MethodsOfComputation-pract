#include <Matrix.hpp>
#include <benchmark/benchmark.h>

static void BM_genfromtxt(benchmark::State &state) {
    for (auto _ : state)
        matrix.genfromtxt<std::string>("./benchmarks/datasets/boston/boston.csv", ',');
}
BENCHMARK(BM_genfromtxt);

BENCHMARK_MAIN();
