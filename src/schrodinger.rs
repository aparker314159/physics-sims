use crate::utils::normalize;


const TOLERANCE: f64 = 1e-5;

// Calculates a solution to the time-independent 
// Schrodinger equation in a symmetric 1D potential.
// Uses Numerov's method, as described in Giannozzi's book
pub struct SchrodingerNumerovSymmetric {
    x_min: f64,
    #[allow(dead_code)]
    x_max: f64, 
    grid_points: usize, // number of points on grid
    pub dx: f64,

}

impl SchrodingerNumerovSymmetric {
    pub fn new(x_max: f64, grid_points: usize) -> Self {
        SchrodingerNumerovSymmetric { 
            x_min: -x_max,
            x_max: x_max,
            grid_points,
            dx: x_max / (grid_points as f64),
        }
    }

    pub fn v(&self, x: f64) -> f64 {
        0.5 * x * x
        //0.5 * ((x / 2.0).powi(4) - 2.0 * (x / 2.0).powi(2) + 1.0)
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
            println!("e = {e}");
            let (r, found_n, d) = self.calculate_with_e(e, n);
            if its > 50 {
                return (r, e);
            }
            println!("{found_n} {d}");
            println!("({e_lower}, {e_upper})");
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
                (e - self.v(x)) * 2.0
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


