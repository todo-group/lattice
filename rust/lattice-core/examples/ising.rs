use lattice_core::Graph;

#[derive(Default)]
struct Observable {
    count: usize,
    sum: f64,
    sum_sq: f64,
}

impl Observable {
    fn add(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.sum_sq += value * value;
    }

    fn mean(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum / self.count as f64
        }
    }

    fn error(&self) -> f64 {
        if self.count < 2 {
            0.0
        } else {
            let n = self.count as f64;
            let mean = self.mean();
            let var = (self.sum_sq / n) - (mean * mean);
            (var.max(0.0) / (n - 1.0)).sqrt()
        }
    }
}

struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_f64(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1);
        let value = self.state >> 11;
        (value as f64) / ((1u64 << 53) as f64)
    }
}

fn main() {
    println!("Metropolis Algorithm for Classical Ferromagnetic Ising Model");

    let beta = 1.0 / 2.269;

    let dim = 2usize;
    let length = 32usize;
    let graph = Graph::simple(dim, length);

    let mut rng = Lcg::new(12345);
    let mut spins = vec![1.0f64; graph.num_sites()];

    let mut energy = Observable::default();
    let mut mag2 = Observable::default();

    let sweeps = 65536usize;
    let therm = sweeps / 8;

    for mcs in 0..(therm + sweeps) {
        for site in 0..graph.num_sites() {
            let mut diff = 0.0;
            for k in 0..graph.num_neighbors(site) {
                diff += 2.0 * spins[site] * spins[graph.neighbor(site, k)];
            }
            if rng.next_f64() < (-beta * diff).exp() {
                spins[site] = -spins[site];
            }
        }

        if mcs > therm {
            let mut ene = 0.0;
            for bond in 0..graph.num_bonds() {
                ene -= spins[graph.source(bond)] * spins[graph.target(bond)];
            }
            energy.add(ene / graph.num_sites() as f64);

            let mag = spins.iter().sum::<f64>() / graph.num_sites() as f64;
            mag2.add(mag * mag);
        }
    }

    println!("dimension = {dim}");
    println!("length = {length}");
    println!("beta = {beta}");
    println!(
        "energy density = {} +/- {}",
        energy.mean(),
        energy.error()
    );
    println!(
        "magnetization density^2 = {} +/- {}",
        mag2.mean(),
        mag2.error()
    );
}
