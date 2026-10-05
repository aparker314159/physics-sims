mod graph;
mod schrodinger;
mod utils;

use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let d = 7.0;
    let a = 0.1;
    schrodinger::run_schrodinger(|x| d * ((-2.0 * a * x).exp() - 2.0 * (-a * x) + 1.0))
}
