use std::error::Error;
use num::complex::Complex;

use cairo::{Context, Format, ImageSurface};
use minifb::{Key, Window, WindowOptions};
use crate::{graph::{self, Graph2d}, utils::{l2product, normalize}};


const TOLERANCE: f64 = 1e-5;
const PRECISION_LEVELS: usize = 20;

pub trait Potential1D : Fn(f64) -> f64 {}
impl <V: Fn(f64)->f64> Potential1D for V {}

// Calculates a solution to the time-independent 
// Schrodinger equation in a symmetric 1D potential.
// Uses Numerov's method, as described in Giannozzi's book

pub struct SchrodingerNumerov<V: Potential1D> {
    x_min: f64,
    #[allow(dead_code)]
    x_max: f64, 
    grid_points: usize, // number of points on grid
    pub dx: f64,
    v: V,
}

impl<V: Potential1D> SchrodingerNumerov<V> {
    pub fn new(v: V, x_max: f64, grid_points: usize) -> Self {
        SchrodingerNumerov::<V> { 
            v,
            x_min: -x_max,
            x_max: x_max,
            grid_points,
            dx: x_max / (grid_points as f64),
        }
    }

    pub fn calculate(&self, n: u64) -> (Vec<(f64, f64)>, f64) {
        /*
        let (r, _, _) = self.calculate_with_e(n as f64 + 0.5, n);
        return r;
        */
        
        let mut e_lower = 0.0;
        let mut e_upper = 50.0;
        let mut e: f64;
        let mut its = 0;

        loop {
            its += 1;
            e = (e_lower + e_upper) / 2.0;
            let (r, found_n, d) = self.calculate_with_e(e, n);
            if its > 50 {
                return (r, e);
            }
            if found_n > n {
                e_upper = e;
                continue;
            } else if found_n < n {
                e_lower = e;
                continue;
            }
            
            if d.abs() < TOLERANCE {
                return (r, e);
            } else if d > 0.0 {
                e_upper = e;
            } else {
                e_lower = e;
            }

            if e_lower == e_upper {
                return (r, e);
            }

        }

    }

    fn calculate_with_e(&self, e: f64, n: u64) -> (Vec<(f64, f64)>, u64, f64) {
        let g: Vec<f64> = (0..=self.grid_points).map(
            |i| {
                let x: f64 = (i as f64) * self.dx;
                (e - (self.v)(x)) * 2.0
            }).collect();

        let f: Vec<f64> = g.clone().into_iter().map(
            |g_i| { 1.0 + g_i * self.dx * self.dx / 12.0 }
        ).collect();

        let mut y_forward: Vec<f64>= vec![];

        let mut meeting_point: usize = 0;
        while g[meeting_point] > 0.0 {
            meeting_point += 1;
            if meeting_point > self.grid_points {
                return (vec![], 10000, 1000.0);
            }
        }

        if n % 2 == 0 {
            let y_0 = self.dx;  // set y_0 to arbitrary value
            y_forward.push(y_0);
            y_forward.push( (12.0 - 10.0 * f[0]) * y_0 / (2.0 * f[1]) );
        } else {
            y_forward.push(0.0);
            y_forward.push(1e-5); // arbitrary
        }

        for i in 2..=meeting_point {
            y_forward.push(
                ((12.0 - 10.0 * f[i-1]) * y_forward[i-1] - f[i-2] * y_forward[i-2]) / f[i]
            );
        }

        /* invert so endpoints are positive */
        if *y_forward.last().expect("") < 0.0 {
            for y_i in &mut y_forward {
                *y_i = -*y_i;
            }
        }

        let mut y_backward: Vec<f64> = vec![];
        y_backward.push(0.0);
        y_backward.push(0.00001); // forces positivity, value is arbitrary
                                  // we're rescaling anyways
        let mut i = self.grid_points - 2;
        let mut j = 2;
        
        while i >= meeting_point {
            y_backward.push(
                ((12.0 - 10.0 * f[i + 1]) * y_backward[j-1] 
                 - f[i + 2] * y_backward[j-2]) 
                / f[i]
            );
            i -= 1;
            j += 1;
        }

        let scale_factor = y_forward.last().expect("") / y_backward.last().expect("");
        for v in &mut y_backward {
            *v *= scale_factor;
        }

        let backward_steps = y_backward.len();

        if backward_steps < 2 {
            let mut y_prev = 1.0;
            let mut sign_changes: u64 = 0;
            for y_i in &y_forward {
                if (*y_i > 0.0 && y_prev < 0.0) || (*y_i < 0.0 && y_prev > 0.0) {
                    sign_changes += 1;
                }
                y_prev = *y_i;
            }
            return (vec![], sign_changes, 1000.0);
        }
        
        let mut y_half = vec![];
        for y_i in &y_forward {
            y_half.push(*y_i);
        }

        y_backward.reverse();
        for y_i in y_backward.iter().skip(1) {
            y_half.push(*y_i);
        }


        let mut y = vec![];
        for y_i in y_half.iter().skip(1).rev() {
            if n % 2 == 0 {
                y.push(*y_i);
            } else {
                y.push(-*y_i);
            }
        }

        for y_i in y_half.iter() {
            y.push(*y_i);
        }
        
        normalize(&mut y, self.dx);

        let new_meeting_point = meeting_point + y_half.len() - 1;
        let discontinuity = (y[new_meeting_point-1] + y[new_meeting_point+1] 
                            - (14.0 - 12.0 * f[meeting_point])*y[new_meeting_point]) 
                            / self.dx;

        let x_vals = (-(self.grid_points as i64)..=(self.grid_points as i64)).map(|i| {(self.x_min + (i as f64)) * self.dx});

        let mut y_prev = y[0];
        let mut sign_changes: u64 = 0;
        for y_i in &y {
            if (*y_i > 0.0 && y_prev < 0.0) || (*y_i < 0.0 && y_prev > 0.0) {
                sign_changes += 1;
            }
            if *y_i != 0.0 {
                y_prev = *y_i;
            }
        }

        return (x_vals.zip(y).collect(), sign_changes, discontinuity);
    }

}


pub fn run_schrodinger<V: Potential1D> (potential: V) -> Result<(), Box<dyn Error>> {
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
    let solver = SchrodingerNumerov::new(&potential,  x_bounds.1, 1000);
    let mut eigenfunctions = vec![];
    for i in 0..=PRECISION_LEVELS {
        eigenfunctions.push(solver.calculate(i as u64));
    }



    let y_bounds: (f64, f64) = (-1.2, 1.2);

    let mut graph = Box::new(Graph2d::new(x_bounds, y_bounds, &mut surface, &context)?);

    let mut start_wavefunction = vec![];

    let mut v_graph = vec![];
    for i in -1000..=1000 {
        let x = i as f64 * solver.dx;
        start_wavefunction.push((x, 1.0/(4.0*x).cosh() + 0.2/(x-1.0).cosh()));
        v_graph.push((x, (potential)(x)));
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
