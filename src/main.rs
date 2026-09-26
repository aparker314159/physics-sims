use std::error::Error;
use num::complex::Complex;

use crate::{graph::Graph2d, schrodinger::SchrodingerNumerovSymmetric, utils::l2product};

mod graph;
mod schrodinger;
mod utils;
use cairo::{Context, Format, ImageSurface};
use minifb::{Key, Window, WindowOptions};


const PRECISION_LEVELS: usize = 20;

fn main() -> Result<(), Box<dyn Error>> {
    let mut window = Window::new(
        "plot",
        graph::WIDTH,
        graph::HEIGHT,
        WindowOptions::default(),
    )?;

    window.set_target_fps(60);
    let mut surface = ImageSurface::create(Format::Rgb24, graph::WIDTH as i32, graph::HEIGHT as i32)?;

    let context = Context::new(&surface)?;
    context.set_antialias(cairo::Antialias::None);

    let x_bounds = (-8.0, 8.0);
    let solver = SchrodingerNumerovSymmetric::new(x_bounds.1, 1000);
    let mut eigenfunctions = vec![];
    for i in 0..=PRECISION_LEVELS {
        eigenfunctions.push(solver.calculate(i as u64));
    }

//    println!("{:?}", data);


    let y_bounds: (f64, f64) = (-1.2, 1.2);

    let mut graph = Box::new(Graph2d::new(x_bounds, y_bounds, &mut surface, &context)?);

    let mut start_wavefunction = vec![];

    let mut v_graph = vec![];
    for i in -1000..=1000 {
        let x = i as f64 * solver.dx;
        start_wavefunction.push((x, 1.0/(4.0*x).cosh() + 0.2/(x-1.0).cosh()));
        v_graph.push((x, solver.v(x)));
    }
    
    let mut c = vec![];
    print!("c = ");
    for i in 0..=PRECISION_LEVELS {
        let ci = l2product(&eigenfunctions[i].0, &start_wavefunction, solver.dx);
        print!("{ci} ");
        c.push(ci);
    }

    let mut t: f64 = 0.0;
    let x_vals: Vec<f64> = eigenfunctions[0].0.iter().map( |(x, _)| {*x} ).collect();
    while window.is_open() && !window.is_key_down(Key::Escape) {

        let mut wavefunction = vec![Complex::new(0.0, 0.0); 2001];

        for i in 0..start_wavefunction.len() {
            for j in 0..=PRECISION_LEVELS {
                wavefunction[i] += c[j] * Complex::new(0.0, -t * eigenfunctions[j].1 ).exp() * eigenfunctions[j].0[i].1
            }
        }


        let re: Vec<f64> = wavefunction.clone()
            .into_iter()
            .map(|v| {v.re})
            .collect();

        /*
        let mut max_diff = 0.0;
        for i in 0..start_wavefunction.len() {
            max_diff = utils::f64_max((re[i] - start_wavefunction[i].1).abs(), max_diff);
        }
        */
        //println!("max = {max_diff}");

        let im: Vec<f64> = wavefunction.clone()
            .into_iter()
            .map(|v| {v.im})
            .collect();

        let prob: Vec<f64> = wavefunction.clone()
            .into_iter()
            .map(|v| {v.norm_sqr()})
            .collect();
        
        let mut graphs = vec![];
        graphs.push(x_vals.clone().into_iter().zip(re).collect());
        graphs.push(x_vals.clone().into_iter().zip(im).collect());
        graphs.push(x_vals.clone().into_iter().zip(prob).collect());
        //graphs.push(v_graph.clone());

        graph.draw_frame(graphs)?;

        window.update_with_buffer(graph.get_pixels(), graph::WIDTH, graph::HEIGHT)?;
        t += 1.0 / 180.0;
    }

    Ok(())
}
