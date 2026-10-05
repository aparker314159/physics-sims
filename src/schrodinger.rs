use std::error::Error;
use num::complex::Complex;

use cairo::{Context, Format, ImageSurface};
use minifb::{Key, Window, WindowOptions};
use crate::{graph::{self, Graph2d}, schrodinger::SchrodingerNumerovErr::{TooHighEnergy, TooLowEnergy}, utils::{count_sign_changes, l2product, normalize, rescale_slice}};


const TOLERANCE: f64 = 1e-5;
const PRECISION_LEVELS: usize = 10;

pub trait Potential1D : Fn(f64) -> f64 {}
impl <V: Fn(f64)->f64> Potential1D for V {}

// Calculates a solution to the time-independent 
// Schrodinger equation in a symmetric 1D potential.
// Uses Numerov's method, as described in Giannozzi's book

pub struct SchrodingerNumerov<V: Potential1D> {
    #[allow(dead_code)]
    x_min: f64,
    #[allow(dead_code)]
    x_max: f64, 
    precision: usize, // number of points on grid
    pub dx: f64,
    v: V,
}

enum SchrodingerNumerovErr {
    TooLowEnergy, 
    TooHighEnergy,
}

impl<V: Potential1D> SchrodingerNumerov<V> {
    pub fn new(v: V, x_max: f64, precision: usize) -> Self {
        SchrodingerNumerov::<V> { 
            v,
            x_min: -x_max,
            x_max: x_max,
            precision, // will have 2*scale + 1 grid points, from -x_max to 0 to x_max.
            dx: x_max / (precision as f64), 
        }
    }

    pub fn calculate(&self, n: u64) -> (Vec<(f64, f64)>, f64) {
        /*
        let (r, _, _) = self.calculate_with_e(n as f64 + 0.5, n);
        return r;
        */
        
        let mut e_lower = 0.0;
        let mut e_upper = 500.0;
        let mut e: f64;
        let mut its = 0;

        loop {
            its += 1;
            e = (e_lower + e_upper) / 2.0;
            let result = self.calculate_with_e(e, n);
            //if let Ok((r, found_n, d)) =  {
            match result {
                Ok((r, found_n, d)) => {
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
                Err(err) => {
                    if its > 50 {
                        panic!("couldn't converge!");
                    } 
                    match err {
                        TooLowEnergy => { e_lower = e },
                        TooHighEnergy => { e_upper = e },
                    }
                }
            }
        }

    }

    fn iterate_numerov_r(y: &mut Vec<f64>, f: &Vec<f64>, start: usize, end: usize) {
        if start <= end {
            for i in start..=end {
                y[i] = ((12.0 - 10.0 * f[i-1]) * y[i-1] 
                       - f[i-2] * y[i-2]) 
                       / f[i];
            }
        }
    }

    fn iterate_numerov_l(y: &mut Vec<f64>, f: &Vec<f64>, start: usize, end: usize) {
        if end <= start {
            let mut i = start;
            while i >= end {
                y[i] = ((12.0 - 10.0 * f[i + 1]) * y[i+1] 
                       - f[i + 2] * y[i+2]) 
                       / f[i];
                i -= 1;
            }
        }
    }

    fn calculate_with_e(&self, e: f64, n: u64) -> Result<(Vec<(f64, f64)>, u64, f64), SchrodingerNumerovErr> {
        let (y_pos, discont_pos)  = self.calculate_half_with_e(e, n, false)?;
        /* 
         * The negative half of the graph can be done the same way by just flipping the sign of the
         * x value fed to the potential, then reversing it at the end.
         */
        let (mut y_neg, _discont_neg) = self.calculate_half_with_e(e, n, true)?;

        /* Glue the two parts together, rescaling the left so the derivative at 0 is continuous.
         * We'll still have two values of discontinuity at the meeting points, but they should
         * share the same sign (I think) */
        let g_0: f64 = e - (self.v)(0.0);
        let g_1: f64 = e - (self.v)(self.dx);
        let g_n1: f64 = e - (self.v)(self.dx * -1.0);
        let f_0: f64 = 1.0 + g_0 * self.dx * self.dx / 12.0;
        let f_1: f64 = 1.0 + g_1 * self.dx * self.dx / 12.0;
        let f_n1: f64 = 1.0 + g_n1 * self.dx * self.dx / 12.0;

        // the expected value of y_neg[1], based on the data from y_pos.
        let expected_y = 
                ((12.0 - 10.0 * f_0) * y_pos[0] 
                 - f_1 * y_pos[1]) 
                / f_n1;

        let rescale_factor = expected_y / y_neg[1];
        rescale_slice(&mut y_neg, rescale_factor, 1, self.precision);
        //println!("rescale factor: {rescale_factor}");

        // Build one large vector with our values.
        let mut y = vec![];
        for y_i in y_neg.into_iter().skip(1).rev() {
            y.push(y_i);
        }
        for y_i in y_pos.into_iter() {
            y.push(y_i);
        }

        normalize(&mut y, self.dx);

        let x: Vec<f64> = (-(self.precision as i64)..=(self.precision as i64))
            .map(|i| i as f64 * self.dx)
            .collect();


        let sign_changes = count_sign_changes(&y);
        return Ok((x.into_iter().zip(y).collect(), sign_changes, discont_pos));
    }

    fn calculate_half_with_e(&self, e: f64, n: u64, going_left: bool) -> Result<(Vec<f64>, f64), SchrodingerNumerovErr> {
        let g: Vec<f64> = (0..=self.precision).map(
            |i| {
                let x: f64 = (i as f64) * self.dx;
                (e - (self.v)(x * if going_left {-1.0} else {1.0})) * 2.0
            }).collect();

        let f: Vec<f64> = g.clone().into_iter().map(
            |g_i| { 1.0 + g_i * self.dx * self.dx / 12.0 }
        ).collect();

        let mut y: Vec<f64>= vec![];
        y.resize(self.precision+1, 0.0);

        let endpoint: usize = self.precision;

        /*
         * FORWARDS FROM ZERO TO MEETING POINT
         */
        let mut meeting_point: usize = 1;
        // note the < here even though g[self.precision] is valid. That's because we index
        // meeting_point+2 later.
        while meeting_point < self.precision-1 && g[meeting_point] > 0.0 {
            meeting_point += 1;
        }

        if meeting_point >= self.precision-1 {
            if (self.v)(self.precision as f64 * self.dx * if going_left {-1.0} else {1.0}) > e {
                return Err(TooLowEnergy);
            } else {
                return Err(TooHighEnergy);
            }
        }


        if n % 2 == 0 {
            let y_0 = self.dx;  // set y_0 to arbitrary value
            y[0] = y_0;
            y[1] = (12.0 - 10.0 * f[0]) * y_0 / (2.0 * f[1]);
        } else {
            y[0] = 0.0;
            y[1] = self.dx;
        }

        Self::iterate_numerov_r(&mut y, &f, 2, meeting_point);
        /*
        for i in 2..=meeting_point_r {
            y_forward_from_0.push(
                ((12.0 - 10.0 * f[i-1]) * y_forward_from_0[i-1] - f[i-2] * y_forward_from_0[i-2]) / f[i]
            );
        }*/
        /* invert so endpoints are positive */
        if y[meeting_point] < 0.0 {
            rescale_slice(&mut y, -1.0, 0, meeting_point);
        }


        /*
         * BACKWARDS FROM RIGHT TO MEETING POINT
         */

        y[endpoint] = 0.0;
        y[endpoint-1] = self.dx; // forces positivity, value is arbitrary
                                   // since we're rescaling anyways
        Self::iterate_numerov_l(&mut y, &f, endpoint - 2, meeting_point+1);
        
        // What does the backwards iteration say about the value at the meeting point?
        let backwards_meeting_point_value: f64 = 
                ((12.0 - 10.0 * f[meeting_point + 1]) * y[meeting_point+1] 
                 - f[meeting_point + 2] * y[meeting_point+2]) 
                / f[meeting_point];

        let scale_factor = y[meeting_point] / backwards_meeting_point_value;
        rescale_slice(&mut y, scale_factor, meeting_point+1, endpoint);


        let discontinuity = (y[meeting_point-1] + y[meeting_point+1] 
                            - (14.0 - 12.0 * f[meeting_point])*y[meeting_point]) 
                            / (self.dx * if going_left { -1.0 } else { 1.0 } );

        return Ok((y, discontinuity));
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
        let solution = solver.calculate(i as u64);
        println!("Energy {i}: {}", solution.1);
        eigenfunctions.push(solution);
    }

    let x_vals: Vec<f64> = eigenfunctions[0].0.iter().map( |(x, _)| {*x} ).collect();
    /*
    let y_bounds: (f64, f64) = (-1.2, 1.2);
    let mut graph = Box::new(Graph2d::new(x_bounds, y_bounds, &mut surface, &context)?);


    while window.is_open() && !window.is_key_down(Key::Escape) {
        let mut graphs = vec![];
        graphs.push(eigenfunctions[0].0.clone());
        graph.draw_frame(graphs)?;
        window.update_with_buffer(graph.get_pixels(), graph::WIDTH, graph::HEIGHT)?;
    }*/

    let y_bounds: (f64, f64) = (-1.2, 1.2);

    let mut graph = Box::new(Graph2d::new(x_bounds, y_bounds, &mut surface, &context)?);

    let mut start_wavefunction = vec![];

    let mut v_graph = vec![];
    for i in -1000..=1000 {
        let x = i as f64 * solver.dx;
        start_wavefunction.push((x, 0.2));
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
