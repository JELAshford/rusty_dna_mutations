use rayon::prelude::*;
use std::time::Instant;

use std::f64::consts::PI;

use ode_solvers::dopri5::*;
use ode_solvers::*;

type State = Vector6<f64>;
type Time = f64;

struct KeplerOrbit {
    mu: f64,
}
impl ode_solvers::System<State> for KeplerOrbit {
    // Equations of motion
    fn system(&self, _t: Time, y: &State, dy: &mut State) {
        let r = (y[0] * y[0] + y[1] * y[1] + y[2] + y[2]).sqrt();
        dy[0] = y[3];
        dy[1] = y[4];
        dy[2] = y[5];
        dy[3] = -self.mu * y[0] / r.powi(3);
        dy[4] = -self.mu * y[1] / r.powi(3);
        dy[5] = -self.mu * y[2] / r.powi(3);
    }
}

fn main() {
    // let mut arr = [0, 1, 2, 3, 4, 5];
    // let timer = Instant::now();
    // arr.par_iter_mut().for_each(|x| *x -= 1);
    // println!("Elased: {:.2?}", &timer.elapsed());

    // Solve example Kepler ode with ode_solvers, see if rayon speedup possible
    // for many configurations
    let config_system = KeplerOrbit { mu: 398600.435436 };
    let a: f64 = 20000.0;
    let period: f64 = 2.0 * PI * (a.powi(3) / config_system.mu).sqrt();
    let y0: State = State::new(
        -5007.248417988539,
        -1444.918140151374,
        3628.534606178356,
        0.717716656891,
        -10.224093784269,
        0.748229399696,
    );

    let start_conditions: Vec<State> = vec![y0; 1_000];

    let timer = Instant::now();
    let _out: Vec<_> = start_conditions
        .par_iter()
        .map(|start_state| {
            let this_system = KeplerOrbit { mu: 398600.435436 };
            let mut stepper = Dopri5::new(
                this_system,
                0.0,
                5.0 * period,
                60.0,
                *start_state,
                1.0e-10,
                1.0e-10,
            );
            stepper.integrate().expect("Integration failed!");
        })
        .collect();
    println!("Elapsed: {:.2?}", &timer.elapsed());
}
