use ode_solvers::dopri5::*;
use ode_solvers::*;

use rand::{rngs::ThreadRng, thread_rng};
use rayon::prelude::*;
use std::collections::HashMap;
use std::time::Instant;

use rand_distr::{Distribution, Normal};

use plotters::prelude::*;
use polars::prelude::*;

type State = Vector4<f64>;
type Time = f64;

struct MutationKinetics {
    kag: f64,
    kat: f64,
    kac: f64,
    kga: f64,
    kgt: f64,
    kgc: f64,
    kta: f64,
    ktg: f64,
    ktc: f64,
    kca: f64,
    kcg: f64,
    kct: f64,
}
impl MutationKinetics {
    fn load_trek_rates() -> HashMap<&'static str, Normal<f64>> {
        include_str!("../data/mutation_rates.txt")
            .lines()
            .map(|line| {
                let parts: Vec<_> = line.split(',').collect();
                (
                    parts[0],
                    Normal::new(
                        parts[1].parse::<f64>().expect("Mean was not a valid float"),
                        parts[2].parse::<f64>().expect("SD was not a valid float"),
                    )
                    .expect("Invalid components"),
                )
            })
            .collect()
    }
    fn from_mutation_rates(
        mutation_rates: &HashMap<&str, Normal<f64>>,
        rng: &mut ThreadRng,
    ) -> MutationKinetics {
        let mut_samples: HashMap<&str, f64> = mutation_rates
            .iter()
            .map(|(name, gen)| (*name, f64::max(gen.sample(rng), 0.0)))
            .collect();
        MutationKinetics {
            kag: mut_samples["AG"],
            kat: mut_samples["AT"],
            kac: mut_samples["AC"],
            kga: mut_samples["GA"],
            kgt: mut_samples["GT"],
            kgc: mut_samples["GC"],
            kta: mut_samples["TA"],
            ktg: mut_samples["TG"],
            ktc: mut_samples["TC"],
            kca: mut_samples["CA"],
            kcg: mut_samples["CG"],
            kct: mut_samples["CT"],
        }
    }
}
impl ode_solvers::System<State> for MutationKinetics {
    fn system(&self, _t: Time, y: &State, dy: &mut State) {
        // y and y0 are ordered AGTC
        dy[0] = self.kca * y[3] + self.kta * y[2] + self.kga * y[1]
            - (self.kac + self.kat + self.kag) * y[0];
        dy[1] = self.kag * y[0] + self.kcg * y[3] + self.ktg * y[2]
            - (self.kga + self.kgt + self.kgc) * y[1];
        dy[2] = self.kat * y[0] + self.kgt * y[1] + self.kct * y[3]
            - (self.kta + self.ktc + self.ktg) * y[2];
        dy[3] = self.kac * y[0] + self.ktc * y[2] + self.kgc * y[1]
            - (self.kca + self.kct + self.kcg) * y[3];
    }
}

fn main() -> anyhow::Result<()> {
    // Analysis structure:
    // - Define Key Parameters
    // - Load Organism-specific starting state
    // - Load Organism/Theory dependant mutation rates
    // - Run paramterised model
    // - Generate plots

    // Define key parameters for runs
    let repeats: usize = 10_000;
    let span: f64 = 5.0;
    let step: f64 = 0.01;
    let organism: &str = "Homo_sapiens";

    // Load and sample base-pair starting states
    let lf1 = LazyFrame::scan_parquet("data/combined_data.parquet", Default::default())?
        .filter(col("Species_name").str().contains(organism) )
        .select([col("A_content"), col("C_content"), col("G_content"), col("T_content")])
        .collect()?;  
    let y0: State = State::new(
        lf1["A_content"].f64()?.to_ndarray()?[0],
        lf1["G_content"].f64()?.to_ndarray()?[0],
        lf1["T_content"].f64()?.to_ndarray()?[0],
        lf1["C_content"].f64()?.to_ndarray()?[0]
    );
    let start_conditions: Vec<State> = vec![y0; repeats];


    // Load mutation rates
    let mut_rates = MutationKinetics::load_trek_rates();


    // Run ODE with sampled rates, timed
    let timer = Instant::now();
    let final_states: Vec<_> = start_conditions
        .par_iter()
        .map(|start_state| {
            let mut this_rng = thread_rng();
            let this_system = MutationKinetics::from_mutation_rates(&mut_rates, &mut this_rng);
            let mut stepper =
                Dopri5::new(this_system, 0.0, span, step, *start_state, 1.0e-10, 1.0e-10);
            stepper.integrate().expect("Integration failed!");
            let ode_out = stepper.y_out().clone();
            ode_out[&ode_out.len() - 1]
        })
        .collect();
    println!("Elapsed: {:.2?}", &timer.elapsed());


    // Extract counts of AGTC from all final states
    let mut a_prop = vec![0.0; repeats];
    let mut g_prop = vec![0.0; repeats];
    let mut t_prop = vec![0.0; repeats];
    let mut c_prop = vec![0.0; repeats];
    let mut at_skew = vec![0.0; repeats];
    let mut gc_skew = vec![0.0; repeats];
    for (ind, row) in final_states.iter().enumerate() {
        a_prop[ind] = row[(0)];
        g_prop[ind] = row[(1)];
        t_prop[ind] = row[(2)];
        c_prop[ind] = row[(3)];
        at_skew[ind] = (row[0] - row[2]) / (row[0] + row[2]);
        gc_skew[ind] = (row[1] - row[3]) / (row[1] + row[3]);
    }


    // Visualise on a scatter plot
    let root = BitMapBackend::new("out/skew_scatter.png", (1024, 1024)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .caption("Skew Distribution", ("sans-serif", 50).into_font())
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-1.1..1.1, -1.1..1.1)?;
    chart.configure_mesh().x_labels(9).y_labels(9).draw()?;
    let mut combo: Vec<(f64, f64)> = vec![(0.0, 0.0); repeats];
    for ind in 0..repeats {
        combo[ind] = (at_skew[ind], gc_skew[ind]);
    }
    chart
        .draw_series(combo.iter().map(|point| Circle::new(*point, 5, RED)))
        .unwrap();
    root.present()?;


    Ok(())
}
