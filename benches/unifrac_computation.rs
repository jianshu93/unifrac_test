use criterion::{
    criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration,
};
use ndarray::Array1;
use ndarray_rand::{
    rand_distr::{Exp1, Uniform},
    RandomExt,
};
use unifrac::compute::{
    parallel_elementwise_sum, parallel_elementwise_sum2, vectorized_elementwise_sum,
};

pub fn bench_rayon_vs_ndarray(c: &mut Criterion) {
    let mut group = c.benchmark_group("UniFrac Function");
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);

    for n_branches in [
        10,
        100,
        1_000,
        10_000,
        100_000,
        1_000_000,
        10_000_000,
        100_000_000,
    ] {
        let p_a = Array1::<f64>::random(n_branches, Uniform::new(0., 1.)).round();
        let p_b = Array1::<f64>::random(n_branches, Uniform::new(0., 1.)).round();
        let brlens = Array1::<f64>::random(n_branches, Exp1);

        group.bench_with_input(
            BenchmarkId::new("External Rayon Parallel", n_branches),
            &(&p_a, &p_b, &brlens),
            |b, (p_a, p_b, brlens)| b.iter(|| parallel_elementwise_sum(p_a, p_b, brlens)),
        );

        group.bench_with_input(
            BenchmarkId::new("Internal Rayon Parallel", n_branches),
            &(&p_a, &p_b, &brlens),
            |b, (p_a, p_b, brlens)| b.iter(|| parallel_elementwise_sum2(p_a, p_b, brlens)),
        );

        group.bench_with_input(
            BenchmarkId::new("NDarray vectorized", n_branches),
            &(&p_a, &p_b, &brlens),
            |b, (p_a, p_b, brlens)| b.iter(|| vectorized_elementwise_sum(p_a, p_b, brlens)),
        );
    }
}

criterion_group!(benches, bench_rayon_vs_ndarray);
criterion_main!(benches);
