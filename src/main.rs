mod graph;
mod schrodinger;
mod utils;

use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    schrodinger::run_schrodinger(|x| 0.5*x*x)
}
