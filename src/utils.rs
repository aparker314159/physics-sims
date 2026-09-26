
pub fn normalize(v: &mut Vec<f64>, dx: f64) {
    /* Integrate v^2 with trapezoid rule */
    let mut norm: f64 = 0.0;

    for i in 0..(v.len()) {
        if i == 0 || i == v.len()-1 {
            norm += 0.5 * v[i] * v[i] * dx;
        } else {
            norm += v[i] * v[i] * dx;
        }
    }
    norm = norm.sqrt();

    for vi in v.into_iter() {
        *vi = *vi / norm;
    }

    println!("norm: {norm}");
}

pub fn l2product(f: &Vec<(f64, f64)>, g: &Vec<(f64, f64)>, dx: f64) -> f64 {
    let mut v: f64 = 0.0;

    for i in 0..(std::cmp::min(f.len(), g.len())) {
        if i == 0 || i == f.len()-1 {
            v += 0.5 * f[i].1 * g[i].1 * dx;
        } else {
            v += f[i].1 * g[i].1 * dx;
        }
    }

    v
}

pub fn f64_max(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}
