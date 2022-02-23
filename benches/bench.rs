use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration,
};
use is_prime_for_primitive_int::IsPrime;
use ring_algorithm::*;
type Z = rug::Integer;

fn make_input(len: usize) -> (Vec<Z>, Vec<Z>) {
    let mut u = Vec::with_capacity(len);
    let mut m = Vec::with_capacity(len);
    for i in 3.. {
        if i.is_prime() {
            let a = rand::random::<u64>() % i;
            u.push(Z::from(a));
            m.push(Z::from(i));
            if m.len() >= len {
                break;
            }
        }
    }
    (u, m)
}

pub fn bench(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("crt");
    group.plot_config(plot_config);
    for i in 5..=15 {
        let size = 1 << i;
        let (u, m) = make_input(size);
        if i <= 11 {
            group.bench_with_input(BenchmarkId::new("crt", size), &size, |b, &_| {
                b.iter(|| chinese_remainder_theorem(black_box(&u), black_box(&m)))
            });
        }
        group.bench_with_input(BenchmarkId::new("fcrt", size), &size, |b, &_| {
            b.iter(|| fast_chinese_remainder_theorem(black_box(&u), black_box(&m)))
        });
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench
}
criterion_main!(benches);
